### This file contains a few useful functions for converting SIMIND output to STIR format.
### It should probably be cleaned up and object-orientedified at some point.

### Author: Sam Porter

import os
import numpy as np
from numbers import Number

from sirf.STIR import ImageData, AcquisitionData
import subprocess
import re
import warnings

def parse_sinogram(template_sinogram):

    template_sinogram.write('tmp.hs')
    values = parse_interfile('tmp.hs')
    os.remove('tmp.hs')
    return values

def parse_interfile(filename):
    """
    Parses STIR interfile file and returns dictionary of key value pairs

    Args:
        filename (string): file name of interfile file

    Returns:
        dict: dictionary of key value pairs
    """
    
    values = {}
    with open(filename, 'r') as file:
        for line in file:
            # Use regex to find lines with ':='
            match = re.search(r'([^;].*?)\s*:=\s*(.*)', line)
            if match:
                key = match.group(1).strip()
                value = match.group(2).strip()
                values[key] = value
    return values

def create_window_file(lower_bounds:list, upper_bounds:list, scatter_orders:list, output_filename:str='input', energy_window=None, lower_ew=None, upper_ew=None):
    """
    Creates a window file for simind simulation

    Args:
        lower_bounds (list): lower bounds of energy windows
        upper_bounds (list): upper bounds of energy windows
        scatter_orders (list): scatter orders of energy windows
        output_filename (str, optional): name of output file. Defaults to 'input'.
        ! energy_window (str, optional): energy window type can be dew or dew. Defaults to None. 
        ! lower_ew (list, optional): lower energy window lower and upper bounds. Defaults to None.
        ! upper_ew (list, optional): upper energy window lower and upper bounds. Defaults to None.
        ! Note that dual and triple energy windows are not yet supported by this wrapper. Please define your own energy windows and work out yourself
    """

    # if path suffix is not. win, add it
    if not output_filename.endswith('.win'):
        output_filename += '.win'

    if isinstance(lower_bounds, Number):
        lower_bounds = [lower_bounds]
    if isinstance(upper_bounds, Number):
        upper_bounds = [upper_bounds]
    if isinstance(scatter_orders, Number):
        scatter_orders = [scatter_orders]

    assert len(lower_bounds) == len(upper_bounds) == len(scatter_orders), "lower_bounds, upper_bounds and scatter_orders must have same length"

    # remove previous window file if present
    if os.path.exists(output_filename):
        os.remove(output_filename)

    with open(output_filename, 'w') as file:
        for i in range(len(lower_bounds)):
            # for some reason simind doesn't like the last line to have a newline character
            file.write(f'{float(lower_bounds[i])},{float(upper_bounds[i])},{int(scatter_orders[i])}\n')
        # simind sometimes doesn't output scatter files unless there's one line with a dedicated scatter order
        # annoyingle this will create one extra total, air, scatter file
        if max(scatter_orders)<1:
            file.write(f'{float(lower_bounds[i])},{float(upper_bounds[i])},1')          
        
def get_sirf_attenuation_from_simind(attn_filename, photopeak_energy = 0.12, attn_type='mu'):
    """Reads attenuation data from simind attenuation file and returns SIRF ImageData object

    Args:
        attn_filename (string): file name of simind attenuation file header
        new_filename (string): file name of new attenuation file header
        attn_type (str, optional): unit of attenuation binary file. Defaults to 'mu'.

    Returns:
        image: SIRF ImageData object containing attenuation data
    """    

    if attn_type == 'mu':
        data_type = np.float32
    elif attn_type == 'rho*1000':
        data_type = np.uint16

    # remove suffix from filename if present
    if attn_filename[-3:] == 'ict' or attn_filename[-3:] == 'hct':
        attn_filename = attn_filename[:-4]

    attn = np.fromfile(attn_filename+".ict", dtype=data_type)
    image = ImageData()

    header_dict = parse_interfile(attn_filename+".hct")
    dim = [int(header_dict['!matrix size [%d]' % i]) for i in range(1,4)][::-1]

    vsize = [float(header_dict['scaling factor (mm/pixel) [%d]' % i]) for i in range(1,3)]
    vsize.append(header_dict['# scaling factor (mm/pixel) [3]'])

    #origin looks like '-128.0000 -128.0000 128.000000'
    origin_string = header_dict['# Image Position First image']
    origin = [float(i) for i in origin_string.split(' ')][::-1]

    image.initialise(dim=tuple(dim), vsize=tuple(vsize), origin=tuple(origin))
    attn = attn.reshape(dim)

    if attn_type == 'mu':
        image.fill(attn)
    elif attn_type == 'rho*100':
        image.fill(attn/1000*photopeak_energy)
    
    return image

def get_sirf_sinogram_from_simind(simind_header_filepath, script_path = ".", circular=True):

    """
    
    ** DEPRECATED METHOD - may still be useful down the line so retained for now **


    Reads sinogram data from simind header file and returns SIRF AcquisitionData object
    Outputs sinogram with counts/MBq/s as units

    Args:
        simind_header_filepath (string): file name of simind header file
        script_path (string, optional): path to simind conversion scripts. Defaults to ".".
        circular (bool, optional): whether to use circular or non-circular conversion script. Defaults to True.
    Returns:
        AcquisitionData: SIRF AcquisitionData object containing sinogram data
    """    

    if circular:
        script = os.path.join(script_path, "convertSIMINDToSTIR.sh")
    else:
        script = os.path.join(script_path, "convertSIMINDToSTIR_noncirc.sh")
    subprocess.run(["sh", os.path.join(script_path, script),\
                     simind_header_filepath])
    
    return AcquisitionData(simind_header_filepath[:-4]+".hs")

def extract_attributes_from_stir_sinogram(sinogram: AcquisitionData):
    
    sinogram.write('tmp.hs')
    attributes = extract_attributes_from_stir_headerfile('tmp.hs')
    #os.remove('tmp.hs')
    return attributes

def extract_attributes_from_stir_headerfile(filename: str):
    attributes = {
        'matrix_sizes': {},
        'scaling_factors': {},
    }

    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().lower()
            if "!imaging modality :=" in line:
                attributes['modality'] = line.split(":=")[1].strip()
            elif "!type of data :=" in line:
                attributes['data_type'] = line.split(":=")[1].strip()
            elif "imagedata byte order :=" in line:
                attributes['byte_order'] = line.split(":=")[1].strip()
            elif "!number format :=" in line:
                attributes['number_format'] = line.split(":=")[1].strip()
            elif "!number of bytes per pixel :=" in line:
                attributes['bytes_per_pixel'] = int(line.split(":=")[1].strip())
            elif "calibration factor:=" in line:
                attributes['calibration_factor'] = float(line.split(":=")[1].strip())
            elif "isotope name:=" in line:
                attributes['isotope_name'] = line.split(":=")[1].strip()
            elif "number of dimensions :=" in line:
                attributes['number_of_dimensions'] = int(line.split(":=")[1].strip())
            elif "!matrix size [" in line:
                axis = line.split('[')[1].split(']')[0]
                attributes['matrix_sizes'][axis] = int(line.split(":=")[1].strip())
            elif "!scaling factor (mm/pixel) [" in line:
                axis = line.split('[')[1].split(']')[0]
                attributes['scaling_factors'][axis] = float(line.split(":=")[1].strip())
            elif "!number of projections :=" in line:
                attributes['number_of_projections'] = int(line.split(":=")[1].strip())
            elif "number of time frames :=" in line:
                attributes['number_of_time_frames'] = int(line.split(":=")[1].strip())
            elif "!image duration (sec)[" in line:
                attributes['image_duration'] = int(line.split(":=")[1].strip())
            elif "!extent of rotation :=" in line:
                attributes['extent_of_rotation'] = float(line.split(":=")[1].strip())
            elif "!direction of rotation :=" in line:
                attributes['direction_of_rotation'] = line.split(":=")[1].strip()
            elif "start angle :=" in line:
                attributes['start_angle'] = float(line.split(":=")[1].strip())
            elif "orbit :=" in line:
                attributes['orbit'] = line.split(":=")[1].strip()
            elif "radius :=" in line or "radii :=" in line:
                tmp = line.split(":=")[1].strip()
                if "{" in tmp and "}" in tmp:  # Handle the list case
                    tmp = tmp.replace("{", "").replace("}", "")
                    values = [float(i) for i in tmp.split(",")]
                    mean_value = np.mean(values)
                    std_value = np.std(values)
                    # Assuming a threshold of 1e-6 for the standard deviation 
                    # to decide if the values are approximately equal
                    if std_value > 1e-6:  
                        attributes['orbit'] = "Non-circular"
                        attributes['radii'] = values
                        warnings.warn("Non-circle orbit detected. There are a number of ways to do this so please handle this case manually using SIMIND Index 42.", UserWarning) 
                    else:
                        attributes['orbit'] = "Circular"
                    attributes['distance_to_detector'] = mean_value
                else:  # Single value case
                    attributes['distance_to_detector'] = float(tmp)
            elif "energy window lower level[1] :=" in line:
                attributes['energy_window_lower'] = float(line.split(":=")[1].strip())
            elif "energy window upper level[1] :=" in line:
                attributes['energy_window_upper'] = float(line.split(":=")[1].strip())
            elif "!name of data file := " in line:
                attributes['data_file'] = line.split(":=")[1].strip()
    return attributes

