import os
import re
import subprocess
import time
from numbers import Number
from simind_attn import *

import numpy as np
from sirf.STIR import AcquisitionData, ImageData

### Should this be more object oriented? ###

def parse_interfile(template_sinogram):

    tmp = template_sinogram.write('tmp.hs')
    values = parse_interfile('tmp.hs')
    os.remove('tmp.hs')
    return values

def parse_interfile(filename):
    """
    Parses interfile file and returns dictionary of key value pairs

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

def create_window_file(lower_bounds:list, upper_bounds:list, scatter_orders:list, output_filename:str='input'):
    """
    Creates a window file for simind simulation

    Args:
        lower_bounds (list): lower bounds of energy windows
        upper_bounds (list): upper bounds of energy windows
        scatter_orders (list): scatter orders of energy windows
        output_filename (str, optional): name of output file. Defaults to 'input'.
    """

    if isinstance(lower_bounds, Number):
        lower_bounds = [lower_bounds]
    if isinstance(upper_bounds, Number):
        upper_bounds = [upper_bounds]
    if isinstance(scatter_orders, Number):
        scatter_orders = [scatter_orders]

    assert len(lower_bounds) == len(upper_bounds) == len(scatter_orders), "lower_bounds, upper_bounds and scatter_orders must have same length"

    # remove previous window file if present
    if os.path.exists(output_filename+".win"):
        os.remove(output_filename+".win")

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

    """Reads sinogram data from simind header file and returns SIRF AcquisitionData object
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
        'modality': None,
        'data_type': None,
        'byte_order': None,
        'number_format': None,
        'bytes_per_pixel': None,
        'calibration_factor': None,
        'isotope_name': None,
        'number_of_dimensions': None,
        'matrix_sizes': {},
        'scaling_factors': {},
        'number_of_projections': None,
        'number_of_time_frames': None,
        'image_duration': None,
        'extent_of_rotation': None,
        'direction_of_rotation': None,
        'start_angle': None,
        'orbit': None,
        'distance_to_detector': None
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
                        raise ValueError("Values are not similar enough for a circular orbit.")
                    attributes['distance_to_detector'] = mean_value
                else:  # Single value case
                    attributes['distance_to_detector'] = float(tmp)
                
    return attributes


def run_simind_simulation(source, mu_map, input_filepath, num_projections, calibration_factor = 1, 
                          photopeak_energy = 0.12, multiplier=1, input_filename = "input", output_filename = "output",
                          output_filepath = ".", source_type="tc99m", collimator = "ma-legp", distance_to_detector = None,
                          start_angle = 360, rotation_angle = 360, direction = "CCW", image_duration_per_projection = 1,
                          source_activity=1):
    """Runs simind simulation - temporary function to inform development of eventual simind_simulator class

    Args:
        source (string or ImageData): source image or file name of source image
        mu_map (string or ImageData): mu map image or file name of mu map image
        input_filepath (string): path to input files
        num_projections (int): number of projections
        calibration_factor (float, optional): calibration factor for source. Defaults to 1.
        photopeak_energy (float, optional): attenuation coefficient of water. Defaults to 0.12.
        multiplier (int, optional): multiplier for number of photons. Defaults to 1.
        input_filename (str, optional): name of input file. Defaults to "input".
        output_filename (str, optional): name of output file. Defaults to "output".
        output_filepath (str, optional): path to output files. Defaults to ".".
        source_type (str, optional): source type. Defaults to "tc99m".
        collimator (str, optional): collimator type. Defaults to "ma-legp".
    """

    if isinstance(source, str):
        source = ImageData(source)
    elif isinstance(source, ImageData):
        pass
    else:
        raise TypeError('source must be a string or SIRF ImageData object')

    if isinstance(mu_map, str):
        mu_map = ImageData(mu_map)
    elif isinstance(mu_map, ImageData):
        pass
    else:
        raise TypeError('mu_map must be a string or SIRF ImageData object')
    
    if source.dimensions() != mu_map.dimensions() or source.voxel_sizes() != mu_map.voxel_sizes():
        raise ValueError('At the moment source and mu_map must have same dimensions and voxel sizes')
    
    source_arr = source.as_array()*calibration_factor
    source_arr.astype(np.float16).tofile('tmp_source.smi')
    
    mu_map_arr = mu_map.as_array()*1000/photopeak_energy
    mu_map_arr.astype(np.uint16).tofile('tmp_density.dmi')

    dim_z, dim_yx, dim_xy = source.dimensions()
    vox_xy, vox_yx, vox_z = source.voxel_sizes()

    #divide voxel size by 10 to get cm
    vox_xy /= 10
    vox_yx /= 10
    vox_z /= 10

    print(vox_xy, vox_yx, vox_z)

    if direction == "CCW":
        rotation_angle = -rotation_angle
        start_angle = -start_angle
    elif direction == "CW":
        pass

    # simind requires a specific number for each direction and amount of rotation
    # STIR can handle any number, but simind can't
    # STIR sinograms can contain floats for start_angle and extent_of_rotation
    # simind sinograms can only contain integers - so we need to convert
    if np.isclose(rotation_angle, -360, rtol=0, atol=1e-1):
        rotation_switch = 0
    elif np.isclose(rotation_angle, -180, rtol=0, atol=1e-1):
        rotation_switch = 1
    elif np.isclose(rotation_angle, 180, rtol=0, atol=1e-1):
        rotation_switch = 2
    elif np.isclose(rotation_angle, 360, rtol=0, atol=1e-1):
        rotation_switch = 3
    else:
        raise ValueError('rotation_angle must be -360, -180, 180 or 360')
    
    start_angle = start_angle + 180
    if start_angle > 360:
        start_angle -= 360

    print(f"rotation_angle: {rotation_angle}")
    print(f"start_angle: {start_angle}")

    if distance_to_detector is None:
        distance_to_detector = vox_xy*dim_xy/2

    if dim_xy != dim_yx or vox_xy != vox_yx:
        raise ValueError('At the moment source and mu map must have square voxels')
    
    for (root, dirs, files) in os.walk(input_filepath):
        for f in files:
            if not os.path.exists(os.path.join(output_filepath,"symlink_"+f)):
                os.symlink(os.path.join(input_filepath,f), os.path.join(output_filepath,"symlink_"+f))
                print("symlink @ " + os.path.join(input_filepath,"symlink_"+f))
            else: print("symlink already exists")

    subprocess.run(["simind", f"symlink_{input_filename}", output_filename,\
                    f"/NN:{multiplier}/FS:tmp_source/FD:tmp_density/FW:symlink_{input_filename}/PX:{vox_xy}\
                    /01:{photopeak_energy} \
                    /TH:{vox_z}/28:{vox_xy}/31:{vox_xy}/fi:{source_type}/cc:{collimator}/ca:1/02:{vox_z*dim_z/2}/05:{vox_z*dim_z/2}\
                    /03:{vox_xy*dim_xy/2}/04:{vox_xy*dim_xy/2}/06:{vox_xy*dim_xy/2}/07:{vox_xy*dim_xy/2}\
                    /14:-1/15:-1/29:{num_projections}/34:{dim_z}/76:{dim_xy}/77:{dim_xy}/78:{dim_xy}\
                    /79:{dim_xy}/81:0/82:0/84:1/30:{rotation_switch}/41:{start_angle}/12:{distance_to_detector}/25:{source_activity}\
                    /tr:04/tr:05/tr:11/tr:14"]) # simulation flags
    
    # remove temporary files
    os.remove('tmp_source.smi')
    os.remove('tmp_density.dmi')

    # remove symlinks
    for (root, dirs, files) in os.walk(output_filepath):
        for f in files:
            if f[:8] == "symlink_":
                try:
                    os.remove(os.path.join(output_filepath,f))
                    print("removed symlink @ " + os.path.join(output_filepath,f))
                except:
                    print("could not remove symlink @ " + os.path.join(output_filepath,f))


### classes ###

import os
import sys


class Converter:
    @staticmethod
    def convert_line(line, dir_switch):

        # lines that require addition of ; prefix
        patterns = [
            ("program", ";"),
            ("patient", ";"),
            ("institution", ";"),
            ("contact", ";"),
            ("ID", ";"),
            ("exam type", ";"),
            ("detector head", ";"),
            ("number of images/energy window", ";"),
            ("time per projection", ";"),
            ("data description", ";"),
            ("total number of images", ";"),
            ("acquisition mode", ";")
        ]
        
        for pattern, prefix in patterns:
            if pattern in line:
                return prefix + line, dir_switch
        
        if "Radius" in line:
            res =  f"Radius := {float(line.split()[2])*10}"
        elif "!number format := short float" in line:
            res = "!number format := float"
        elif "image duration" in line:
            parts = line.split()
            res =  f"number of time frames := 1\nimage duration (sec) [1] := {parts[4]}"
        elif ";energy window lower level" in line:
            res =  f"energy window lower level[1] := {line.split()[5]}"
        elif ";energy window upper level" in line:
            res =  f"energy window upper level[1] := {line.split()[5]}"
        elif "CCW" in line:
            res = line
            dir_switch = -1
        elif "start angle" in line:
            res =  f"start angle := {dir_switch*float(line.split()[3]) + 180}"
        else:
            res = line
        
        return res, dir_switch

    @staticmethod
    def convert(filename, return_object=True):
        dir_switch = 1 # -1 for CCW, 1 for CW
        if not filename.endswith(".h00"):
            print(f"USAGE: {sys.argv[0]} filename.h00")
            sys.exit(1)

        stirfilename = filename.replace(".h00", ".hs")
        with open(filename, 'r') as f_in, open(stirfilename, 'w') as f_out:
            for line in f_in:
                write_line, dir_switch = Converter.convert_line(line.strip(), dir_switch)
                f_out.write(write_line+"\n")
        
        print(f"Output in {stirfilename}")

        if return_object:
            return AcquisitionData(stirfilename)

    @staticmethod
    def fix_values_based_on_reference(reference_file, file_to_adjust, threshold, output_adjusted_file=None):
        """
        Read two files and align values in file_to_adjust based on reference_file if they are within a threshold.
        
        Parameters:
        - reference_file: Path to the reference file.
        - file_to_adjust: Path to the file which values might be adjusted.
        - threshold: The threshold for value differences.
        - output_adjusted_file: Path to save the aligned output.
        """

        if output_adjusted_file is None:
            output_adjusted_file = file_to_adjust[:-4] + "_adjusted.h00"
        
        with open(reference_file, 'r') as ref_file:
            reference_lines = ref_file.readlines()

        with open(file_to_adjust, 'r') as to_adjust_file:
            adjust_lines = to_adjust_file.readlines()

        for index, line_to_adjust in enumerate(adjust_lines):
            if ":=" in line_to_adjust:
                # Extract the potential numeric values after ':='
                key_to_adjust = line_to_adjust.split(":=")[0].strip()
                value_to_adjust_part = line_to_adjust.split(":=")[1].strip()

                # Find the line in reference file
                for line_ref in reference_lines:
                    if key_to_adjust + " :=" in line_ref:
                        value_ref_part = line_ref.split(":=")[1].strip()

                        try:
                            value_ref = float(value_ref_part)
                            value_to_adjust = float(value_to_adjust_part)

                            # If the difference between the values is within the threshold
                            if abs(value_ref - value_to_adjust) <= threshold:
                                adjust_lines[index] = key_to_adjust + " := " + str(value_ref) + '\n'
                                break  # Once we've found a match in the reference, we can break out of the loop

                        except ValueError:  # Non-numeric value encountered
                            pass  # Keep the original line from file_to_adjust

        with open(output_adjusted_file, 'w') as out:
            out.writelines(adjust_lines)

        return AcquisitionData(output_adjusted_file)
    
    @staticmethod
    def replace_values_based_on_reference(reference_file, file_to_adjust, output_adjusted_file=None):
        """
        Read two files and align values in file_to_adjust based on reference_file if they are within a threshold.
        
        Parameters:
        - reference_file: Path to the reference file.
        - file_to_adjust: Path to the file which values might be adjusted.
        - threshold: The threshold for value differences.
        - output_adjusted_file: Path to save the aligned output.
        """

        if output_adjusted_file is None:
            output_adjusted_file = file_to_adjust[:-4] + "_adjusted.h00"
        
        with open(reference_file, 'r') as ref_file:
            reference_lines = ref_file.readlines()

        with open(file_to_adjust, 'r') as to_adjust_file:
            adjust_lines = to_adjust_file.readlines()

        for index, line_to_adjust in enumerate(adjust_lines):
            if ":=" in line_to_adjust:
                # Extract the potential numeric values after ':='
                key_to_adjust = line_to_adjust.split(":=")[0].strip()
                value_to_adjust_part = line_to_adjust.split(":=")[1].strip()

                # Find the line in reference file
                for line_ref in reference_lines:
                    if key_to_adjust + " :=" in line_ref:
                        value_ref_part = line_ref.split(":=")[1].strip()

                        try:
                            value_ref = float(value_ref_part)
                            value_to_adjust = float(value_to_adjust_part)
                            adjust_lines[index] = key_to_adjust + " := " + str(value_ref) + '\n'

                        except ValueError:  # Non-numeric value encountered
                            pass  # Keep the original line from file_to_adjust

        with open(output_adjusted_file, 'w') as out:
            out.writelines(adjust_lines)

        return AcquisitionData(output_adjusted_file)
    
    @staticmethod
    def replace_sinogram_values_based_on_reference(reference_sinogram, sinogram_to_adjust):

        # write tmp files
        reference_sinogram.write('tmp_ref.hs')
        sinogram_to_adjust.write('tmp_adjust.hs')

        # replace values
        result = Converter.replace_values_based_on_reference('tmp_ref.hs', 'tmp_adjust.hs')

        # remove tmp files
        os.remove('tmp_ref.hs')
        os.remove('tmp_adjust.hs')
        os.remove('tmp_ref.s')
        os.remove('tmp_adjust.s')

        return result

    @staticmethod
    def convert_sinogram_parameter(sinogram, parameter, value):
        filename = "tmp.hs"
        sinogram.write(filename)
        # Read the file into memory
        with open(filename, 'r') as f:
            lines = f.readlines()

        # Modify the appropriate line based on the parameter
        with open(filename, 'w') as f:
            for line in lines:
                # Check if the current line starts with the given parameter (ignoring case)
                line_stripped = line.strip().lower()
                param_with_colon = parameter.lower() + " :="
                if line_stripped.startswith(param_with_colon) or line_stripped.startswith('!' + param_with_colon):
                    if line.startswith('!'):
                        line = f"!{parameter} := {value}\n"
                    else:
                        line = f"{parameter} := {value}\n"
                f.write(line)
        sinogram = AcquisitionData(filename)
        os.remove(filename)
        return sinogram


import os
import subprocess

import numpy as np


class SimindSimulator:

    """ Class to run simind simulations
    Currently only for circular orbit, but should be trivial to add non-circular orbit
    """    
    
    def __init__(self, input_filepath, output_filepath="."):
        """ Initialises SimindSimulator object

        Args:
            input_filepath (str): path to simind input files
            output_filepath (str, optional): _description_. Defaults to ".".
        """        
        self.input_filepath = input_filepath
        self.output_filepath = output_filepath

        self.source = None # source image - should be SIRF ImageData object or file name
        self.mu_map = None # mu map image - should be SIRF ImageData object or file name
        self.num_projections = None # number of projections
        self.calibration_factor = 1 # calibration factor for output sinogram
        self.photopeak_energy = 140 # photopeak energy in keV - defaults to Tc-99m
        self.multiplier = 0.01 # multiplier for number of photons simulated
        self.input_filename = "input" # input filename - defaults to input. Must be the same for all input files at the moment
        self.output_filename = "output" # output filename - defaults to output. Is the same for all output files at the moment
        self.source_type = "tc99m" # source type
        self.collimator = "ma-legp" # collimator - defaults to Mediso AnyScan LE GP
        self.distance_to_detector = None
        self.start_angle = 360 # start angle for camera rotation
        self.rotation_angle = 360 # number of degrees to rotate camera
        self.direction = "CCW" # direction of rotation
        self.image_duration_per_projection = 1 # image duration per projection in seconds
        self.source_activity = 1 # source activity in MBq
        self.model_attenuation = True # whether to model attenuation in simind
        self.fix_rounded_output_values = True # whether to fix rounded output values in simind
        self.kev_per_channel = 2 # keV per channel in simind - keV starts at 0
        self.num_energy_spectra_channels = 200 # number of energy spectra channels in simind
        self.template_sinogram = None
        self.output = None
        self.custom_variables = {}
        self.custom_flags = {}

        self.attribute_map = {
            'SourceMap': 'source',
            'MuMap': 'mu_map',
            'input_filepath': 'input_filepath',
            'NumberOfProjections': 'num_projections',
            'CalibrationFactor': 'calibration_factor',
            'PhotopeakEnergy': 'photopeak_energy',
            'PhotonMultiplier': 'multiplier',
            'input_filename': 'input_filename',
            'output_filename': 'output_filename',
            'output_filepath': 'output_filepath',
            'SourceType': 'source_type',
            'Collimator': 'collimator',
            'DistanceToDetector': 'distance_to_detector',
            'StartAngle': 'start_angle',
            'RotationAngle': 'rotation_angle',
            'DirectionOfRotation': 'direction',
            'ImageDurationPerProjection': 'image_duration_per_projection',
            'SourceActivity': 'source_activity',
            'ModelAttenuation': 'model_attenuation',
            'FixRoundedOutputValues': 'fix_rounded_output_values',
            'keVPerChannel': 'kev_per_channel',
            'NumberOfEnergySpectraChannels': 'num_energy_spectra_channels'
        }

        self.files_converted = False # flag to check if output files have been converted to stir

    def set_attributes(self, **kwargs):
        """Set attributes based on keyword arguments"""
        for key, value in kwargs.items():
            self.set_input(key, value)

    def set_attributes(self, attributes):
        """Set attributes based on dictionary of attributes"""
        for key, value in attributes.items():
            self.set_input(key, value)
    
    def get_attributes(self):
        """Returns dictionary of attributes"""
        return {key: getattr(self, value) for key, value in self.attribute_map.items()}
    
    def print_attributes(self):
        """Prints dictionary of attributes with line breaks"""
        for key, value in self.attribute_map.items():
            print(f"{key}: {getattr(self, value)}")

    def check_set_up(self):
        """Check that all required attributes are set before running the simulation."""
        
        required_attributes = self.attribute_map.values()
        
        for attribute in required_attributes:
            if not hasattr(self, attribute):
                raise ValueError(f"The attribute '{attribute}' must be set before running the simulation.")
        
        print("All required attributes are set and ready for simulation.")
        return True

    def check_relevant_files(self):
        """Check that all relevant files exist in the input_filepath."""
        
        assert os.path.exists(os.path.join(self.input_filepath, self.input_filename + ".smc")), f"input_filename.smc does not exist in {self.input_filepath}"
        assert os.path.exists(os.path.join(self.input_filepath, self.input_filename + ".win")), f"input_filename.win does not exist in {self.input_filepath}"

        print(f"All relevant files exist in {self.input_filepath}.")

    def clear_output_filepath(self):
        """Clears all files in output_filepath"""

        print(f"Clearing all files in {self.output_filepath} that could cause trouble...")
        print("This includes:")
        print(" - .h00, .hs, .a00, .hct, .ict, .bis, .res")
        print("if this is not what you want, please cancel this process now.")
        # wait 5 seconds and count down
        for i in range(5,0,-1):
            print(i)
            time.sleep(1)

        for (root, dirs, files) in os.walk(self.output_filepath):
            del dirs[:] # don't recurse into subdirectories
            for f in files:
                if f[-4:] == ".h00" or f[-3:] == ".hs" or f[-4:] == ".a00" or f[-4:] == ".hct" or f[-4:] == ".ict" or f[-4:] == ".bis" or f[-4:] == ".res" and self.output_filename in f:
                    os.remove(os.path.join(root, f))
                    print(f"removed {f}")

    def set_windows(self, lower_bounds, upper_bounds, scatter_orders):
        """Sets energy windows for simind simulation"""
        create_window_file(lower_bounds, upper_bounds, scatter_orders, 
                           output_filename = os.path.join(self.input_filepath, 
                                                          self.input_filename + ".win"))

    def check_source_density_match(self):
        """Checks that source and mu map have same dimensions and voxel sizes"""
        assert self.source.voxel_sizes() == self.mu_map.voxel_sizes(), "Source and mu map must have same voxel sizes"
        assert self.source.dimensions() == self.mu_map.dimensions(), "Source and mu map must have same dimensions"

    def check_square_pixels_and_image(self, image):
        """ Checks that image has square pixels and same x,y dimensions"""

        assert image.voxel_sizes()[2] == image.voxel_sizes()[1], "Image must have square pixels"
        assert image.dimensions()[1] == image.dimensions()[2], "Image must have same x,y dimensions"

    def print_identifiers(self):
        """Print all possible string identifiers for class attributes."""
        for identifier in self.attribute_map:
            print(identifier)

    def set_input(self, identifier, value):
        """Set an attribute based on a string identifier."""

        if identifier not in self.attribute_map:
            raise ValueError(f"'{identifier}' is not a valid identifier.\nValid identifiers can be found by calling print_identifiers()")
        
        setattr(self, self.attribute_map[identifier], value)

    def get_input(self, identifier):
        """Get an attribute based on a string identifier."""

        if identifier not in self.attribute_map:
            raise ValueError(f"'{identifier}' is not a valid identifier.\nValid identifiers can be found by calling print_identifiers()")
        
        return getattr(self, self.attribute_map[identifier])

    def set_source(self, source):
        if isinstance(source, str):
            self.source = ImageData(source)
        elif isinstance(source, ImageData):
            self.source = source
        else:
            raise TypeError('source must be a string or SIRF ImageData object')

    def set_mu_map(self, mu_map):
        if isinstance(mu_map, str):
            self.mu_map = ImageData(mu_map)
        elif isinstance(mu_map, ImageData):
            self.mu_map = mu_map
        else:
            raise TypeError('mu_map must be a string or SIRF ImageData object')
        
    def set_model_attenuation(self, model_attenuation: bool):
        self.model_attenuation = model_attenuation

    def set_num_projections(self, num_projections):
        self.num_projections = num_projections

    def set_calibration_factor(self, calibration_factor):
        self.calibration_factor = calibration_factor

    def set_photopeak_energy(self, photopeak_energy:float):
        self.photopeak_energy = photopeak_energy

    def set_multiplier(self, multiplier: Number):
        self.multiplier = multiplier

    def set_input_filename(self, input_filename: str):
        self.input_filename = input_filename

    def set_output_filename(self, output_filename: str):
        self.output_filename = output_filename

    def set_source_type(self, source_type: str):
        self.source_type = source_type

    def set_collimator(self, collimator: str):
        self.collimator = collimator

    def set_distance_to_detector(self, distance_to_detector: Number):
        self.distance_to_detector = distance_to_detector

    def set_start_angle(self, start_angle: Number):
        self.start_angle = start_angle

    def set_rotation_angle(self, rotation_angle: Number):
        self.rotation_angle = int(np.round(rotation_angle, 0))

    def set_direction(self, direction: str):
        direction = direction.upper().strip()
        if direction not in ["CW", "CCW"]:
            raise ValueError("direction must be 'CW' or 'CCW'")
        self.direction = direction

    def set_image_duration_per_projection(self, image_duration_per_projection: Number):
        self.image_duration_per_projection = image_duration_per_projection

    def set_source_activity(self, source_activity: Number):
        """ Multiplies the values of detector bins by this factor."""
        self.source_activity = source_activity
        
    def set_custom_variable(self, variable_code, value):
        """Sets a custom variable for the Simind simulation."""
        self.custom_variables[variable_code] = value
        
    def set_custom_flag(self, flag_code, value:bool):
        """Sets a custom flag for the Simind simulation."""
        if value:
            value_str = "tr"
        else:
            value_str = "fa"
        self.custom_flags[value_str] = flag_code

    def set_template_sinogram(self, template_sinogram):
        """ Sets the template sinogram for the simulation. """
        print("Warning: This will overwrite any other settings with those found in the template sinogram.")
        if isinstance(template_sinogram, str):
            attribute_dict = extract_attributes_from_stir_headerfile(template_sinogram)
        elif isinstance(template_sinogram, AcquisitionData):
            attribute_dict = extract_attributes_from_stir_sinogram(template_sinogram)
        else:
            raise TypeError('template_sinogram must be a string or SIRF AcquisitionData object')
        self.num_projections = attribute_dict['number_of_projections']
        self.start_angle = attribute_dict['start_angle']
        self.rotation_angle = int(np.round(attribute_dict['extent_of_rotation'], 0))
        self.direction = attribute_dict['direction_of_rotation']
        self.distance_to_detector = attribute_dict['distance_to_detector']/10

        self.template_sinogram = template_sinogram.clone()

    def run_simulation(self):
        """
        Runs simind simulation - method for SimindSimulator class
        Currently only supports square pixels. 
        """

        self.files_converted = False # flag to check if output files have been converted to stir
        self.clear_output_filepath()

        self.check_set_up()
        self.check_relevant_files()

        if isinstance(self.source, str):
            self.source = ImageData(self.source)
        elif not isinstance(self.source, ImageData):
            raise TypeError('source must be a string or SIRF ImageData object')

        if isinstance(self.mu_map, str):
            self.mu_map = ImageData(self.mu_map)
        elif not isinstance(self.mu_map, ImageData):
            raise TypeError('mu_map must be a string or SIRF ImageData object')
        
        if self.source.dimensions() != self.mu_map.dimensions() or self.source.voxel_sizes() != self.mu_map.voxel_sizes():
            raise ValueError('At the moment source and mu_map must have same dimensions and voxel sizes')
        
        self.check_source_density_match() # currently only supports pixels of same size in source and mu_map
        self.check_square_pixels_and_image(self.source)
        self.check_square_pixels_and_image(self.mu_map)

        # save temporary binary files
        source_arr = self.source.as_array() * self.calibration_factor
        source_arr.astype(np.float16).tofile('tmp_source.smi')

        if self.model_attenuation:
            mu_map_arr = self.mu_map.as_array()
            mu_map_arr = attenuation_to_density(mu_map_arr, self.photopeak_energy)*1000
            import matplotlib.pyplot as plt
            mu_map_arr.astype(np.uint16).tofile('tmp_density.dmi')
            att_switch = "tr"
        else:
            att_switch = "fa"

        # get dimensions and voxel sizes
        dim_z, _, dim_xy = self.source.dimensions()
        vox_xy, _, vox_z = self.source.voxel_sizes()

        # Divide voxel size by 10 to get cm
        vox_xy /= 10
        vox_z /= 10

        if self.direction.lower() == "ccw":
            # simind needs angles to be negative for CCW rotation
            rotation_angle = -self.rotation_angle
            start_angle = -self.start_angle
        elif self.direction.lower() == "cw":
            rotation_angle = self.rotation_angle
            start_angle = self.start_angle
        else:
            raise ValueError("direction must be 'CW' or 'CCW'")

        # simind requires a specific number for each direction and amount of rotation
        # STIR can handle any number, but simind can't
        # STIR sinograms can contain floats for start_angle and extent_of_rotation
        # simind sinograms can only contain integers - so we need to convert
        if np.isclose(rotation_angle, -360, rtol=0, atol=1e-1):
            rotation_switch = 0
        elif np.isclose(rotation_angle, -180, rtol=0, atol=1e-1):
            rotation_switch = 1
        elif np.isclose(rotation_angle, 360, rtol=0, atol=1e-1):
            rotation_switch = 2
        elif np.isclose(rotation_angle, 180, rtol=0, atol=1e-1):
            rotation_switch = 3
        else:
            raise ValueError('rotation_angle must be -360, -180, 180 or 360')
        
        # STIR and simind have different conventions for start angle
        start_angle = self.start_angle + 180
        if start_angle > 360:
            start_angle -= 360

        # if distance_to_detector is not set, set it to half the source map size
        if self.distance_to_detector is None:
            self.distance_to_detector = vox_xy * dim_xy / 2
        
        for (root, dirs, files) in os.walk(self.input_filepath):
            del dirs[:] # don't recurse into subdirectories
            for f in files:
                if not os.path.exists(os.path.join(self.output_filepath, "symlink_" + f)) and self.input_filename in f:
                    os.symlink(os.path.join(self.input_filepath, f), os.path.join(self.output_filepath, "symlink_" + f))
                else:
                    print("symlink already exists")

        command = [
            "simind", f"symlink_{self.input_filename}", self.output_filename, # input and output filenames
            f"/NN:{self.multiplier}" # photon multiplier (accuracy versus time)
            f"/FS:tmp_source/FD:tmp_density/FW:symlink_{self.input_filename}" # source, density and energy window filenames
            f"/14:-1/15:-1" # set to use voxelised phantoms
            f"/fi:{self.source_type}" # SIMIND source type (e.g Lu177 / Tc99m)
            f"/cc:{self.collimator}" # define collimator
            "/84:1" # scattwin routine
            f"/01:{self.photopeak_energy}" # photopeak energy
            f"/PX:{vox_xy}/TH:{vox_z}" # pixel size and slice thickness
            f"/28:{vox_xy}/31:{vox_xy}" # pixel size in x and y (I think unnecessary due to above)
            f"/02:{vox_z * dim_z / 2}/05:{vox_z * dim_z / 2}/03:{vox_xy * dim_xy / 2}\
                /04:{vox_xy * dim_xy / 2}/06:{vox_xy * dim_xy / 2}/07:{vox_xy * dim_xy / 2}" # source and mumap sizes
            f"/29:{self.num_projections}/34:{dim_z}/76:{dim_xy}/77:{dim_xy}/78:{dim_xy}/79:{dim_xy}" # number of projections and dimensions of simulated data, source and phantom
            f"/81:0/82:0" # use above (78,79) for source & density matrix sizes
            f"/30:{rotation_switch}/41:{start_angle}/12:{self.distance_to_detector}" # rotation switch, start angle and distance to detector
            f"/80:{self.num_energy_spectra_channels}/27:{self.kev_per_channel}" # energy spectra channels and keV per channel
            f"/25:{self.image_duration_per_projection*self.source_activity}" # source activity * image duration
            "/tr:04" # include collimator
            "/tr:05"  # simulate SPECT (as opposed to planar)
            f"/{att_switch}:11" # whether to include attenuation
            "/tr:14" # set to write interfile header
        ]
        
        # Add custom variables
        for variable_code, value in self.custom_variables.items():
            command.append(f"/{variable_code}:{value}")
            
        # Add custom flags
        for flag_code, value in self.custom_flags.items():
            command.append(f"/{value}:{flag_code}")
            
        subprocess.run(command)

        # remove temporary files
        os.remove('tmp_source.smi')
        os.remove('tmp_density.dmi')

        # remove symlinks
        for (root, dirs, files) in os.walk("."):
            del dirs[:] # don't recurse into subdirectories
            for f in files:
                if f.startswith("symlink_"):
                    try:
                        os.remove(os.path.join(self.output_filepath, f))
                        print("removed symlink @ " + os.path.join(self.output_filepath, f))
                    except:
                        print("could not remove symlink @ " + os.path.join(self.output_filepath, f))


    def get_output(self):
        """Get output files from simind simulation"""
        # convert to .hs files
        
        if self.output is not None:
            return self.output
        
        converter = Converter()
        if not self.files_converted:
            # find all files with output directory ending in .h00 
            files = [f for f in os.listdir(self.output_filepath) if f.endswith('.h00')]
            for f in files:
                converter.convert(os.path.join(self.output_filepath, f))
            self.files_converted = True
        # find all files with output directory ending in .hs
        files = [f for f in os.listdir(self.output_filepath) if f.endswith('.hs')]
        # sort converted files by number (w{n}.hs) and then tot_w{n}.hs, sca_w{n}.hs, air_w{n}.hs
        # files look like {output_filename}_{tot/sca/air}_w{n}.hs
        files.sort(key=lambda f: int(f.split('_')[-1].split('.')[0][1:]))
        # return dictionary of files with keys tot, sca, air and w{n}
        output = {}
        for f in files[:-1]:
            f_split = f.split('_')
            scat_type = f_split[-2]
            window = f_split[-1].split('.')[0]
            output[scat_type + "_" + window] = AcquisitionData(os.path.join(self.output_filepath, f))
            # simind can round these values so let's make sure they are correct
            if self.template_sinogram is not None:
                output[scat_type + "_" + window] = converter.replace_sinogram_values_based_on_reference(self.template_sinogram, output[scat_type + "_" + window])
            else:
                converter.convert_sinogram_parameter(output[scat_type + "_" + window], "scaling factor (mm/pixel) [1]", self.source.voxel_sizes()[1])
                converter.convert_sinogram_parameter(output[scat_type + "_" + window], "scaling factor (mm/pixel) [2]", self.source.voxel_sizes()[2])
                converter.convert_sinogram_parameter(output[scat_type + "_" + window], "Radius", self.distance_to_detector*10)
        
        self.output = output
        
        return output
    
    def get_output_total(self, window=1):
        """Get total output file from simind simulation"""
        outputs = self.get_output()
        return outputs['tot_w' + str(window)]
    
    def get_output_scatter(self, window=1):
        """Get scatter output file from simind simulation"""
        outputs = self.get_output()
        return outputs['sca_w' + str(window)]
    
    def get_output_air(self, window=1):
        """Get air output file from simind simulation"""
        outputs = self.get_output()
        return outputs['air_w' + str(window)]
