import os
import re
import subprocess
from pathlib import Path
import time
from numbers import Number
from .simind_attn import attenuation_to_density

import numpy as np
from sirf.STIR import AcquisitionData, ImageData

### Should this be more object oriented? ###

def parse_sinogram(template_sinogram):

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
                        raise ValueError("Values are not similar enough for a circular orbit.")
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
            ("acquisition mode", ";"),

        ]
        
        for pattern, prefix in patterns:
            if pattern in line:
                return prefix + line, dir_switch
        
        if "Radius" in line:
            res =  f"Radius := {float(line.split()[3])*10}"
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
        elif "!name of data file" in line:
            file = Path(line.split()[5])
            # only filename plus extension
            res = f"!name of data file := {file.stem+file.suffix}"
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

        # close files
        f_in.close()
        f_out.close()

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

    """ 
    Class to run simind simulations with SIRF data
    Currently only for circular orbit, but should be trivial to add non-circular orbit
    """    
    
    def __init__(self, template_smc_file_path, output_dir, output_prefix='output', source=None, mu_map=None, template_sinogram=None, photon_multiplier=1):
        """ Initialises SimindSimulator object

        Args:
            input_filepath (str): path to simind input files
            output_dir (str): path to output directory
        """

        # Ensure the input template file exists
        self.check_files_exist([template_smc_file_path])

        # Convert the template file path to a Path object
        self.template_smc_file_path = Path(template_smc_file_path)
        self.input_dir = self.template_smc_file_path.parent

        # Define the path for the output SMC file
        self.smc_file_path = self.input_dir / "simind.smc"

        # Ensure the output directory exists, create it if not
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Define the path for the output file with the given prefix
        self.output_filepath = self.output_dir / output_prefix

        self.config = SimulationConfig(template_smc_file_path)
        self.runtime_switches = RuntimeSwitches()

        self.source = source
        self.mu_map = mu_map
        self.template_sinogram = template_sinogram
        
        self.output = None

        if source is not None:
            self.set_source(source)
        if mu_map is not None:
            self.set_mu_map(mu_map)
        if template_sinogram is not None:
            self.set_template_sinogram(template_sinogram)

        self.runtime_switches.set_switch('NN', photon_multiplier)

        # set appropriate indices and flags for voxelised phantom
        self.config.set_flag(5, True) # SPECT study
        self.config.set_value(15, -1) # source type
        self.config.set_value(14, -1) # phantom type
        self.config.set_flag(14, True) # write to interfile header
        
        self.window_set = False

    def check_files_exist(self, filepaths):
        """Check that file exists"""
        for f in filepaths:
            assert os.path.exists(f), f"{f} does not exist"

    @staticmethod   
    def find_parent_directory(file_path):
        # Get the absolute path of the file
        absolute_path = os.path.abspath(file_path)
        # Get the directory name of the file path
        directory_name = os.path.dirname(absolute_path)
        # Get the parent directory
        parent_directory = os.path.dirname(directory_name)
        return parent_directory

    # output filename is scattwin.win in the same dir as the template smc file
    def set_windows(self, lower_bounds, upper_bounds, scatter_orders):
        output_filename = os.path.join(self.input_dir, "scattwin.win")
        """Sets energy windows for simind simulation"""
        create_window_file(lower_bounds, upper_bounds, scatter_orders, 
                           output_filename)
        
        self.window_set = True

    def check_images_match(self, image0, image1):
        """Checks that source and mu map have same dimensions and voxel sizes"""
        assert image0.voxel_sizes() == image1.voxel_sizes(), "Source and mu map must have same voxel sizes"
        assert image0.dimensions() == image1.dimensions(), "Source and mu map must have same dimensions"

    def check_square_pixels_and_image(self, image):
        """ Checks that image has square pixels and same x,y dimensions"""
        assert image.voxel_sizes()[2] == image.voxel_sizes()[1], "Image must have square pixels"
        assert image.dimensions()[1] == image.dimensions()[2], "Image must have same x,y dimensions" 

    def add_index(self, index, value):
        """Add an index value to the simulation"""
        self.config.set_value(index, value)
        
    def add_flag(self, flag, value:bool):
        """Add a flag to the simulation"""
        self.config.set_flag(flag, value)

    def add_text_variable(self, variable, value):
        """Add a text variable to the simulation"""
        self.config.set_text_variable(variable, value)

    def add_data_file(self, data_file):
        """Add a data file to the simulation"""
        self.config.set_data_file(data_file)

    def add_comment(self, comment):
        """Add a comment to the simulation"""
        self.config.set_comment(comment)

    def add_runtime_switch(self, switch, value):
        """Add a runtime switch to the simulation"""
        self.runtime_switches.set_switch(switch, value)

    def set_source(self, source):
        if isinstance(source, str):
            self.source = ImageData(source)
        elif isinstance(source, ImageData):
            self.source = source
        else:
            raise TypeError('source must be a string or SIRF ImageData object')
        
        # get dimensions and voxel sizes
        dim_z, dim_xy, vox_xy, vox_z = self.get_dimensions_and_voxel_sizes(self.source)
        # Divide voxel size by 10 to get cm
        vox_xy /= 10
        vox_z /= 10
        
        self.add_index(2, vox_z*dim_z/2)
        self.add_index(3, vox_xy*dim_xy/2)
        self.add_index(4, vox_xy*dim_xy/2)

        self.add_index(28, vox_xy)
        self.add_index(76, dim_xy)
        self.add_index(77, dim_xy)

    def set_mu_map(self, mu_map):
        if isinstance(mu_map, str):
            self.mu_map = ImageData(mu_map)
        elif isinstance(mu_map, ImageData):
            self.mu_map = mu_map
        else:
            raise TypeError('mu_map must be a string or SIRF ImageData object')
        
        # get dimensions and voxel sizes
        dim_z, dim_xy, vox_xy, vox_z = self.get_dimensions_and_voxel_sizes(self.mu_map)

        # Divide voxel size by 10 to get cm
        vox_xy /= 10
        vox_z /= 10
        
        self.add_index(5, vox_z*dim_z/2)
        self.add_index(6, vox_xy*dim_xy/2)
        self.add_index(7, vox_xy*dim_xy/2)

        self.add_index(31, vox_xy) # pixel size density images
        self.add_index(33, 1) # first image density images
        self.add_index(34, dim_z) # number density images
        self.add_index(78, dim_xy) # matrix size density map i
        self.add_index(79, dim_xy) # matrix size density map j

        self.runtime_switches.set_switch('PX', vox_xy)
        self.runtime_switches.set_switch('TH', vox_z)
        
    def set_rotation_in_stir_geometry(self, rotation_angle:Number, start_angle:Number, direction:str,):

        if direction.lower() == "ccw":
            # simind needs angles to be negative for CCW rotation
            rotation_angle = -rotation_angle
            start_angle = -start_angle
        elif direction.lower() == "cw":
            rotation_angle = rotation_angle
            start_angle = start_angle
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
        
        # start angles in simind and STIR are opposite
        start_angle+=180
        if start_angle > 360:
            start_angle -= 360
        
        return rotation_switch, start_angle

    def set_template_sinogram(self, template_sinogram):
        """ Sets the template sinogram for the simulation. """
        print("Warning: This will overwrite any other settings with those found in the template sinogram.")
        if isinstance(template_sinogram, str):
            attribute_dict = extract_attributes_from_stir_headerfile(template_sinogram)
        elif isinstance(template_sinogram, AcquisitionData):
            attribute_dict = extract_attributes_from_stir_sinogram(template_sinogram)
        else:
            raise TypeError('template_sinogram must be a string or SIRF AcquisitionData object')
        
        self.add_index(29, attribute_dict['number_of_projections'])
        self.add_index(41, attribute_dict['start_angle'])
        self.add_index(12, attribute_dict['distance_to_detector']/10) # convert to cm
        rotation_switch, start_angle = self.set_rotation_in_stir_geometry(attribute_dict['extent_of_rotation'], 
                                                                          attribute_dict['start_angle'],
                                                                          attribute_dict['direction_of_rotation'],)
        self.add_index(30, rotation_switch)
        self.add_index(41, start_angle)

        self.template_sinogram = template_sinogram.clone()

    @staticmethod
    def update_linux_path_strings(path):
        """Updates path strings for Linux by converting slashes to backslashes"""

        return path.as_posix().replace("/", "\\")

    @staticmethod
    def reset_linux_path_strings(path):
        """Resets path strings to POSIX style and converts them back to Path objects"""

        return Path(path.replace("\\", "/"))

    def get_dimensions_and_voxel_sizes(self, image):
        """Get dimensions and voxel sizes of image"""
        dim_z, _, dim_xy = image.dimensions()
        vox_xy, _, vox_z = image.voxel_sizes()
        return dim_z, dim_xy, vox_xy, vox_z

    def run_simulation(self):
        """
        Runs simind simulation - method for SimindSimulator class
        Currently only supports square pixels. 
        """

        self.files_converted = False # flag to check if output files have been converted to stir
        
        if self.window_set is False:
            raise ValueError("Energy windows must be set before running simulation\nUse set_windows method")
        
        self.check_images_match(self.source, self.mu_map)
        self.check_square_pixels_and_image(self.source)
        self.check_square_pixels_and_image(self.mu_map)

        
        if self.config.get_flag(11):
            mu_map_arr = self.mu_map.as_array()
            mu_map_arr = attenuation_to_density(mu_map_arr, self.config.get_value('photon_energy'), self.input_dir)*1000
        else:
            mu_map_arr = np.zeros(self.mu_map.shape)
            
        mu_map_arr.astype(np.uint16).tofile('tmp_density.dmi')
        self.config.set_data_file(11, "tmp_density")
            
        self.source.as_array().astype(np.float16).tofile('tmp_source.smi')
        self.config.set_data_file(12, "tmp_source")

        # write smc file
        self.config.save_file(self.output_filepath)

        # if linux os update path strings
        if os.name == 'posix':
            print("Updating path strings for linux")
            output_filepath = self.update_linux_path_strings(self.output_filepath)
        else:
            output_filepath = self.output_filepath

        command = ["simind", output_filepath, output_filepath]

        # add switches
        switches = ""
        for key, value in self.runtime_switches.switches.items():
            switches+=(f'/{key}:{str(value)}')
        command.append(switches)
            
        print(f"Running simind with command: {' '.join(command)}")     
            
        subprocess.run(command)

        # remove temporary files
        os.remove('tmp_source.smi')
        os.remove('tmp_density.dmi') 
        
        # check if output files have been put in the output directory
        if len(os.listdir(self.output_dir)) == 0:
            print("No output files found in output directory. SIMIND isn't very good at this\nManually moving files. Sorry if this moves files you don't want moved")
            for f in os.listdir(self.input_dir):
                if f.endswith('.h00') or f.endswith('.a00'):
                    os.rename(os.path.join(self.input_dir, f), os.path.join(self.output_dir, f))

    def get_output(self):
        """Get output files from simind simulation"""
        # convert to .hs files
        
        if self.output is not None and len(self.output) > 0:
            return self.output
        converter = Converter()
        output_strings = ["air", "sca", "tot", "pri"]
        if not self.files_converted:
            # find all files with output directory ending in .h00 
            files = [f for f in os.listdir(self.output_dir) if f.endswith('.h00') and any(s in f for s in output_strings)]
            for f in files:
                converter.convert(os.path.join(self.output_dir, f))
                print(f"Converted {f}")
            self.files_converted = True
        # if output dir is not empty, convert files
        # find all files with output directory ending in .hs
        files = [f for f in os.listdir(self.output_dir) if f.endswith('.hs') and any(s in f for s in output_strings)]
        print(files)
        # sort converted files by number (w{n}.hs) and then tot_w{n}.hs, sca_w{n}.hs, air_w{n}.hs
        # files look like {output_filename}_{tot/sca/air}_w{n}.hs
        files.sort(key=lambda f: int(f.split('_')[-1].split('.')[0][1:]))
        # return dictionary of files with keys tot, sca, air and w{n}
        output = {}
        for f in files[:-1]:
            f_split = f.split('_')
            scat_type = f_split[-2]
            window = f_split[-1].split('.')[0]
            print(os.path.join(self.output_dir, f))
            output[scat_type + "_" + window] = AcquisitionData(os.path.join(self.output_dir, f))
            # simind can round these values so let's make sure they are correct
            if self.template_sinogram is not None:
                output[scat_type + "_" + window] = converter.replace_sinogram_values_based_on_reference(self.template_sinogram, output[scat_type + "_" + window])
            else:
                converter.convert_sinogram_parameter(output[scat_type + "_" + window], "scaling factor (mm/pixel) [1]", self.source.voxel_sizes()[1])
                converter.convert_sinogram_parameter(output[scat_type + "_" + window], "scaling factor (mm/pixel) [2]", self.source.voxel_sizes()[2])
                converter.convert_sinogram_parameter(output[scat_type + "_" + window], "Radius", float(self.config.get_value("distance_to_detector"))*10)
        
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
    

import re

class SimulationConfig:
    def __init__(self, filepath):
        self.filepath = filepath
        self.index_dict = {
            1: "photon_energy", 2: "source_half_length", 3: "source_half_width", 4: "source_half_height",
            5: "phantom_half_length", 6: "phantom_half_width", 7: "phantom_half_height", 8: "crystal_half_length_radius",
            9: "crystal_thickness", 10: "crystal_half_width", 11: "backscattering_material_thickness",
            12: "height_to_detector_surface", 13: "cover_thickness", 14: "phantom_type", 15: "source_type",
            16: "shift_source_x", 17: "shift_source_y", 18: "shift_source_z", 19: "photon_direction",
            20: "upper_window_threshold", 21: "lower_window_threshold", 22: "energy_resolution",
            23: "intrinsic_resolution", 24: "emitted_photons_per_decay", 25: "source_activity",
            26: "number_photon_histories", 27: "kev_per_channel", 28: "pixel_size_simulated_image",
            29: "spect_no_projections", 30: "spect_rotation", 31: "pixel_size_density_images",
            32: "orientation_density_images", 33: "first_image_density_images", 34: "number_density_images",
            35: "density_limit_border", 36: "shift_density_images_y", 37: "shift_density_images_z",
            38: "step_size_photon_path_simulation", 39: "shift_density_images_x", 40: "density_threshold_soft_bone",
            41: "spect_starting_angle", 42: "spect_orbital_rotation_fraction", 43: "camera_offset_x",
            44: "camera_offset_y", 45: "code_definitions_zubal_phantom", 46: "hole_size_x", 47: "hole_size_y",
            48: "distance_between_holes_x", 49: "distance_between_holes_y", 50: "shift_center_hole_x",
            51: "shift_center_hole_y", 52: "collimator_thickness", 53: "collimator_routine", 54: "hole_shape",
            55: "type", 56: "distance_collimator_detector", 57: "unused_parameter_1", 58: "unused_parameter_2",
            59: "random_collimator_movement", 60: "unused_parameter_3", 76: "matrix_size_image_i",
            77: "matrix_size_image_j", 78: "matrix_size_density_map_i", 79: "matrix_size_source_map_i",
            80: "energy_spectra_channels", 81: "matrix_size_density_map_j", 82: "matrix_size_source_map_j",
            83: "cutoff_energy_terminate_photon_history", 84: "scoring_routine", 85: "csv_file_content",
            91: "voltage", 92: "mobility_life_electrons", 93: "mobility_life_holes", 94: "contact_pad_size",
            95: "anode_element_pitch", 96: "exponential_decay_constant_tau", 97: "components_hecht_formula",
            98: "energy_resolution_model", 99: "cloud_mobility", 100: "detector_array_size_i",
            101: "detector_array_size_j"
        }
        self.flag_dict = {
            1: "write_results_to_screen", 2: "write_images_to_files", 3: "write_pulse_height_distribution_to_file",
            4: "include_collimator", 5: "simulate_spect_study", 6: "include_characteristic_xray_emissions",
            7: "include_backscattering_material", 8: "use_random_seed_value", 9: "currently_not_in_use",
            10: "include_interactions_in_cover", 11: "include_interactions_in_phantom",
            12: "include_energy_resolution_in_crystal", 13: "include_forced_interactions_in_crystal",
            14: "write_interfile_header_files", 15: "save_aligned_phantom_images"
        }
        
        self.data_file_dict = {
            7: "phantom_soft_tissue", 8: "phantom_bone", 9: "cover_material", 10: "crystal_material",
            11: "image_file_phantom", 12: "image_file_source", 13: "backscatter_material", 14: "energy_resolution_file",
        }
        
        self.dict = {"index": self.index_dict, "flag": self.flag_dict, "data_file": self.data_file_dict}
        
        self.data = None
        self.flags = None
        self.text_variables = {}
        self.data_files = {}
        self.comment = None
        self.read_file()


    def read_file(self):
        with open(self.filepath, 'r') as file:
            lines = file.readlines()
            self.comment = lines[1].strip()

            # Parsing Basic Change data
            data_lines = lines[3:27]
            data_string = ' '.join(data_lines).replace("\n", "")
            self.data = [float(val) for val in re.findall(r'-?\d+\.\d+E[+-]\d+', data_string)]

            # Parsing Simulation flags
            self.flags = lines[28].strip().replace(' ', '')

            # Parsing Text Variables
            text_variables_start = 29
            text_variables_count = int(lines[text_variables_start].split()[0])
            text_variables_lines = lines[text_variables_start + 1:text_variables_start + 1 + text_variables_count]
            self.text_variables = {i+1: text_variables_lines[i].strip() for i in range(text_variables_count)}

            # Parsing Data files
            data_files_start = 38
            data_files_count = int(lines[data_files_start].split()[0])
            data_files_lines = lines[data_files_start + 1:data_files_start + 1 + data_files_count]
            self.data_files = {i+7: data_files_lines[i].strip() for i in range(data_files_count)}
            
    def print_config(self):
        print(f"Comment: {self.comment}")
        print("Basic Change data:")
        for key, val in self.index_dict.items():
            print(f"index {key}: {val}: {self.data[key - 1]}")
        print("Simulation flags:")
        for key, val in self.flag_dict.items():
            print(f"flag {key}: {val}: {self.flags[key - 1]}")
        print("Text Variables:")
        for key, val in self.text_variables.items():
            print(f"{key}: {val}")
        print("Data Files:")
        for key, val in self.data_files.items():
            print(f"{key}: {val}")

    def _get_value_by_index(self, index):
        return self.data[index - 1]
    
    def _get_value_by_description(self, description):
        for key, val in self.index_dict.items():
            if val == description:
                return self.data[key - 1]
            
    def get_value(self, index):
        if isinstance(index, int) and index in self.index_dict:
            return self._get_value_by_index(index)
        elif isinstance(index, str) and index in self.index_dict.values():
            return self._get_value_by_description(index)
        else:
            raise ValueError("index must be a valid integer or string")

    def _set_value_by_index(self, index, value):
        self.data[index - 1] = value

    def _set_value_by_description(self, description, value):
        for key, val in self.index_dict.items():
            if val == description:
                self.data[key - 1] = value

    def set_value(self, index, value):
        if isinstance(index, int) and index in self.index_dict:
            self._set_value_by_index(index, value)
        elif isinstance(index, str) and index in self.index_dict.values():
            self._set_value_by_description(index, value)
        else:
            raise ValueError("index must be an integer or string")

    def _get_flag_by_index(self, index):
        return self.flags[index - 1] == 'T'
    
    def _get_flag_by_description(self, description):
        for key, val in self.flag_dict.items():
            if val == description:
                return self.flags[key - 1] == 'T'
            
    def get_flag(self, index):
        if isinstance(index, int) and index in self.flag_dict:
            return self._get_flag_by_index(index)
        elif isinstance(index, str) and index in self.flag_dict.values():
            return self._get_flag_by_description(index)
        else:
            raise ValueError("index must be a valid integer or string")

    def _set_flag_by_index(self, index, value):
        flags = list(self.flags)
        flags[index - 1] = 'T' if value else 'F'
        self.flags = ''.join(flags)

    def _set_flag_by_description(self, description, value):
        for key, val in self.flag_dict.items():
            if val == description:
                flags = list(self.flags)
                flags[key - 1] = 'T' if value else 'F'
                self.flags = ''.join(flags)

    def set_flag(self, index, value):
        if isinstance(index, int) and index in self.flag_dict:
            self._set_flag_by_index(index, value)
        elif isinstance(index, str) and index in self.flag_dict.values():
            self._set_flag_by_description(index, value)
        else:
            raise ValueError("index must be an integer or string")

    def get_text_variable(self, var_index):
        return self.text_variables.get(var_index)

    def set_text_variable(self, var_index, value):
        if var_index in self.text_variables:
            self.text_variables[var_index] = value

    def _get_data_file_by_index(self, index):
        return self.data_files.get(index)
    
    def _get_data_file_by_description(self, description):
        for key, val in self.data_files.items():
            if val == description:
                return self.data_files[key]
            
    def get_data_file(self, index):
        if isinstance(index, int) and index in self.data_files:
            return self._get_data_file_by_index(index)
        elif isinstance(index, str) and index in self.data_files.values():
            return self._get_data_file_by_description(index)
        else:
            raise ValueError("index must be a valid integer or string")

    def _set_data_file_by_index(self, index, value):
        if index in self.data_files:
            self.data_files[index] = value
         
    def _set_data_file_by_description(self, description, value):
        for key, val in self.data_files.items():
            if val == description:
                self.data_files[key] = value
                
    def set_data_file(self, index, value):
        if isinstance(index, int) and index in self.data_file_dict:
            self._set_data_file_by_index(index, value)
        elif isinstance(index, str) and index in self.data_file_dict.values():
            self._set_data_file_by_description(index, value)
        else:
            raise ValueError(f"index must be an integer or string") 
            
    def get_comment(self):
        return self.comment

    def set_comment(self, comment):
        self.comment = comment

    def save_file(self, filepath):
        # Ensure filepath is a Path object
        filepath = Path(filepath)
        
        # Check if the file has the correct suffix, add it if missing
        if filepath.suffix != '.smc':
            filepath = filepath.with_suffix('.smc')
        
        # Proceed with saving the file
        # Your file-saving logic here, using the updated `filepath`

        with open(filepath, 'w') as file:
            comment = self.comment + " "*(70 - len(self.comment))
            file.write(f"SMCV2\n{comment}\n")
            file.write("   120  # Basic Change data\n")

            for i in range(0, len(self.data), 5):
                line = ''
                for val in self.data[i:i+5]:
                    # Format the value in scientific notation with 5 decimal places
                    formatted_val = f"{val:.5E}"   
                    if val != 0:
                        # Split the formatted value into its components: sign, digit, and exponent
                        sign = '-' if val < 0 else ' '
                        parts = formatted_val.split('E')
                        digits = parts[0].replace('-', '')
                        # remove final 0
                        digits = digits[:-1]
                        exponent = int(parts[1])

                        # Ensure it starts with '0' after the sign
                        if '.' in digits:
                            digits = digits.replace('.', '')
                        
                        # Since we've moved the decimal place one position to the right, increment the exponent
                        new_exponent = exponent + 1

                        # Reconstruct the formatted value
                        formatted_val = f"{sign}0.{digits}E{new_exponent:+03d}"
                    else:
                        # If the value is 0, we don't need to format it
                        formatted_val = f" {val:.5E}"
                    # Add the formatted value to the line with padding to ensure consistent spacing
                    line += f"{formatted_val}"

                # Write the formatted line to the file
                file.write(f"{line}\n")

            file.write(f"    30  # Simulation flags\n{self.flags}\n")
            file.write(f"     {len(self.text_variables)}  # Text Variables\n")
            for i in range(1, len(self.text_variables) + 1):
                file.write(f"{self.text_variables[i]}\n")
            file.write(f"    {len(self.data_files)} # Data files\n")
            for i in range(7, 7 + len(self.data_files)):
                # needs to be 60 long including spaces
                data_file = self.data_files[i] + " "*(60 - len(self.data_files[i]))
                file.write(f"{data_file}\n")

    # methods for printing what each index and flag corresponds to
    def print_index_dict(self):
        for key, value in self.index_dict.items():
            print(f"{key}: {value}")

    def print_flag_dict(self):
        for key, value in self.flag_dict.items():
            print(f"{key}: {value}")

    def print_index(self, index):
        print(f"{index}: {self.index_dict[index]}")

    def print_flag(self, flag):
        print(f"{flag}: {self.flag_dict[flag]}")
        
    @property
    def combined_dict(self):
        combined_dict = {}
        for sub_dict in self.dict.values():
            combined_dict.update(sub_dict)
        return combined_dict

class RuntimeSwitches:
    def __init__(self):
        self.standard_switch_dict = {
            'CC': 'Collimator code',
            'DF': 'Density file segment',
            'ES': 'Energy offset',
            'FE': 'Energy resolution file',
            'FZ': 'Zubal file',
            'FI': 'Input file',
            'FD': 'Density map base name',
            'FS': 'Source map base name',
            'I2': 'Image files stored as 16-bit integer matrices',
            'IN': 'Change simind.ini value',
            'LO': 'Photon histories before printout',
            'LF': 'Linear sampling of polar angle for photon direction',
            'MP': 'MPI parallel run',
            'OR': 'Change orientation of the density map',
            'PR': 'Start simulation at projection number',
            'PU': 'Shift of the source in pixel units',
            'QF': 'Quit simulation if earlier result file exists',
            'RR': 'Random number generator seed',
            'SC': 'Maximum number of scatter orders',
            'SF': 'Segment for source map',
            'TS': 'Time shift for interfile header',
            'UA': 'Set density equal to data buffer or 1.0',
            'WB': 'Whole-body simulation of anterior and posterior views',
            'Xn': 'Change cross sections'
        }

        self.image_based_switch_dict = {
            'PX': 'Pixel size of the source maps',
            'DI': 'General direction of the source map',
            'TH': 'Slice thickness for the images',
            'SB': 'Start block when reading source maps',
            '1S': 'Position of the first image to be used',
            'NN': 'Multiplier for scaling the number of counts',
            'IF': 'Input tumour file',
        }

        self.myocardiac_switch_dict = {
            'A1': 'Shift of the heart in the xy-direction',
            'A2': 'Shift of the heart in the yz-direction',
            'A3': 'Shift of the heart in the zx-direction',
            'L1': 'Location of defect',
            'L2': 'Angular size of the defect',
            'L3': 'Start from Base',
            'L4': 'Extent of defect in axis direction',
            'L5': 'Transgression in %',
            'L6': 'Activity ratio in defect',
            'M1': 'Thickness of the myocardial wall',
            'M2': 'Thickness of the plastic wall',
            'M3': 'Total length of the chamber',
            'M4': 'Total diameter of the chamber',
        }

        self.multiple_spheres_switch_dict = {
            'C1': 'Number of spheres',
            'C2': 'Radius of spheres',
            'C3': 'Activity of spheres',
            'C4': 'Shift of spheres in the x-direction',
            'C5': 'Shift of spheres in the y-direction',
            'C6': 'Shift of spheres in the z-direction',
        }
        self.switches = {}

        self.switch_dict = {"Standard": self.standard_switch_dict,
                            "Image-based": self.image_based_switch_dict,
                            "Myocardiac": self.myocardiac_switch_dict,
                            "Multiple spheres": self.multiple_spheres_switch_dict}
        
    @property
    def combined_switch_dict(self):
        combined_dict = {}
        for sub_dict in self.switch_dict.values():
            combined_dict.update(sub_dict)
        return combined_dict

    def _set_switch_by_switch(self, switch, value):
        if switch in self.combined_switch_dict:
            self.switches[switch] = value
        else:
            raise ValueError(f"Switch {switch} is not recognised.")
        
    def _set_switch_by_name(self, name, value):
        for switch, description in self.combined_switch_dict.items():
            if description == name:
                self.switches[switch] = value
                return
        raise ValueError(f"Switch {name} is not recognised.")
    
    def set_switch(self, identifier, value):
        if identifier in self.combined_switch_dict.values():
            self._set_switch_by_name(identifier, value)
        elif identifier in self.combined_switch_dict.keys():
            self._set_switch_by_switch(identifier, value)
        else:
            raise ValueError(f"Switch {identifier} is not recognised.")

    def print_switches(self):
        for switch, value in self.switches.items():
            description = self.combined_switch_dict[switch]
            print(f"{switch} ({description}): {value}")

    def print_available_switches(self):
        for switch_dict in self.switch_dict.values():
            print(f"Switches for {switch_dict}:")
            for switch, description in switch_dict.items():
                print(f"{switch}: {description}")