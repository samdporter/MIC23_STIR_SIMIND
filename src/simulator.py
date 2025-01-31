### This file contains classed for SIMIND simulations with SIRF data.
### It reads in a template SMC file and runs a simulation with the given source and mu map.
### The output is then converted to STIR format.

### Author: Sam Porter

import os
import subprocess
from pathlib import Path
import numpy as np
from numbers import Number
from sirf.STIR import ImageData, AcquisitionData
from .converter import Converter

from .simind_attn import attenuation_to_density
from .config import SimulationConfig, RuntimeSwitches
from .functions import extract_attributes_from_stir_headerfile, extract_attributes_from_stir_sinogram, create_window_file


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
        self.template_smc_file_path = self.template_smc_file_path.resolve()
        self.input_dir = self.template_smc_file_path.parent

        # Define the path for the output SMC file
        self.smc_file_path = self.input_dir / "simind.smc"
        self.smc_file_path = self.smc_file_path.resolve()

        # Ensure the output directory exists, create it if not
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir = self.output_dir.resolve()

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
        output_filename = os.path.join(self.output_dir, "scattwin.win")
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
        """ Sets the template sinogram for the simulation. 
        
        This will overwrite any other settings with those found in the template sinogram.
        Settings include:
        - number of projections
        - distance to detector
        - rotation direction
        - start angle (converted to simind geometry)
        
        """
        print("Warning: This will overwrite any other settings with those found in the template sinogram.")
        if isinstance(template_sinogram, str):
            attribute_dict = extract_attributes_from_stir_headerfile(template_sinogram)
        elif isinstance(template_sinogram, AcquisitionData):
            attribute_dict = extract_attributes_from_stir_sinogram(template_sinogram)
        else:
            raise TypeError('template_sinogram must be a string or SIRF AcquisitionData object')
    
        self.add_index(29, attribute_dict['number_of_projections'])
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

        self.output = None

        self.files_converted = False # flag to check if output files have been converted to stir
        
        if self.window_set is False:
            raise ValueError("Energy windows must be set before running simulation\nUse set_windows method")
        
        self.check_images_match(self.source, self.mu_map)
        self.check_square_pixels_and_image(self.source)
        self.check_square_pixels_and_image(self.mu_map)
        
        cwd = os.getcwd()
        os.chdir(self.output_dir)

        if self.config.get_flag(11):
            mu_map_arr = self.mu_map.as_array()
            mu_map_arr = attenuation_to_density(mu_map_arr, self.config.get_value('photon_energy'), self.input_dir)*1000
        else:
            mu_map_arr = np.zeros(self.mu_map.shape)
            
        mu_map_arr = mu_map_arr.astype(np.uint16)
        mu_map_arr.tofile('tmp_density.dmi')
        self.config.set_data_file(11, "tmp_density")
            
            
        # make dynamic range of source 0-100
        source_arr = self.source.as_array()
        source_arr/=self.source.max()
        source_arr*=100
        source_arr = np.round(source_arr).astype(np.uint16)
        source_arr.tofile('tmp_source.smi')
        self.config.set_data_file(12, "tmp_source")

        # write smc file
        self.config.save_file(self.output_filepath)

        command = ["simind", self.output_filepath.name, self.output_filepath.name]

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
                    
        os.chdir(cwd)

    def get_output(self):
        """Get output files from simind simulation"""
        # convert to .hs files
        
        if self.output is not None and len(self.output) > 0:
            return self.output
        converter = Converter()
        output_strings = ["_air_w", "_sca_w", "_tot_w", "_pri_w"]
        if not self.files_converted:
            # find all files with output directory ending in .h00 
            h00_files = [f for f in os.listdir(self.output_dir) if f.endswith('.h00') and any(s in f for s in output_strings)]
            for f in h00_files:
                converter.convert(os.path.join(self.output_dir, f))
                print(f"Converted {f}")
            self.files_converted = True
        # if output dir is not empty, convert files
        # find all files with output directory ending in .hs
        hs_files = [f for f in os.listdir(self.output_dir) if f.endswith('.hs') and any(s in f for s in output_strings)]
        # sort converted files by number (w{n}.hs) and then tot_w{n}.hs, sca_w{n}.hs, air_w{n}.hs
        # files look like {output_filename}_{tot/sca/air}_w{n}.hs
        hs_files.sort(key=lambda f: int(f.split('_')[-1].split('.')[0][1:]))
        # return dictionary of files with keys tot, sca, air and w{n}
        output = {}
        for f in hs_files:
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

        # remove h00 files
        for f in h00_files:
            os.remove(os.path.join(self.output_dir, f))
        
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