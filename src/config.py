### This file contains a wrapper to access, edit and save simulation configuration files for the Simind Monte Carlo simulation software.
### It can work as a standalone to make the .smc files accessible and editable in a more user-friendly way.
### Or you can use it with the Simulator class to run the simulation with SIRF in python.

### Author: Sam Porter

import re
from pathlib import Path

class SimulationConfig:
    """
    SimulationConfig Class

    This class is designed to parse, manipulate, and save simulation configuration files. It provides easy access
    to configuration parameters, including index-based data, simulation flags, text variables, and associated data files.

    Attributes:
        filepath (str): Path to the simulation configuration file.
        index_dict (dict): Dictionary mapping indices to parameter names for basic change data.
        flag_dict (dict): Dictionary mapping indices to simulation flags.
        data_file_dict (dict): Dictionary mapping indices to data file descriptions.
        data (list): List of basic change data values.
        flags (str): String representing simulation flags as 'T' (True) or 'F' (False).
        text_variables (dict): Dictionary of text variables.
        data_files (dict): Dictionary of associated data files.
        comment (str): Comment section from the configuration file.
    """

    def __init__(self, filepath):
        """
        Initialize the SimulationConfig instance.

        Args:
            filepath (str): Path to the simulation configuration file.
        """
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
            11: "image_file_phantom", 12: "image_file_source", 13: "backscatter_material", 14: "energy_resolution_file"
        }
        self.data = None
        self.flags = None
        self.text_variables = {}
        self.data_files = {}
        self.comment = None
        self.read_file()

    def read_file(self):
        """
        Parse the simulation configuration file and populate attributes.
        """
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
        """
        Print the configuration details, including comments, basic change data, flags, text variables, and data files.
        """
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

    def get_value(self, index):
        """
        Get the value of a parameter by its index or description.

        Args:
            index (int or str): Parameter index or description.

        Returns:
            float: Parameter value.
        """
        if isinstance(index, int) and index in self.index_dict:
            return self.data[index - 1]
        elif isinstance(index, str) and index in self.index_dict.values():
            for key, val in self.index_dict.items():
                if val == index:
                    return self.data[key - 1]
        else:
            raise ValueError("index must be a valid integer or string")

    def set_value(self, index, value):
        """
        Set the value of a parameter by its index or description.

        Args:
            index (int or str): Parameter index or description.
            value (float): New value for the parameter.
        """
        if isinstance(index, int) and index in self.index_dict:
            self.data[index - 1] = value
        elif isinstance(index, str) and index in self.index_dict.values():
            for key, val in self.index_dict.items():
                if val == index:
                    self.data[key - 1] = value
        else:
            raise ValueError("index must be an integer or string")

    def get_flag(self, index):
        """
        Get the value of a simulation flag by its index or description.

        Args:
            index (int or str): Flag index or description.

        Returns:
            bool: True if the flag is set, False otherwise.
        """
        if isinstance(index, int) and index in self.flag_dict:
            return self.flags[index - 1] == 'T'
        elif isinstance(index, str) and index in self.flag_dict.values():
            for key, val in self.flag_dict.items():
                if val == index:
                    return self.flags[key - 1] == 'T'
        else:
            raise ValueError("index must be a valid integer or string")

    def set_flag(self, index, value):
        """
        Set the value of a simulation flag by its index or description.

        Args:
            index (int or str): Flag index or description.
            value (bool): True to set the flag, False to clear it.
        """
        if isinstance(index, int) and index in self.flag_dict:
            flags = list(self.flags)
            flags[index - 1] = 'T' if value else 'F'
            self.flags = ''.join(flags)
        elif isinstance(index, str) and index in self.flag_dict.values():
            for key, val in self.flag_dict.items():
                if val == index:
                    flags = list(self.flags)
                    flags[key - 1] = 'T' if value else 'F'
                    self.flags = ''.join(flags)
        else:
            raise ValueError("index must be an integer or string")
        
    def set_data_file(self, index, filepath):
        """
        Set the path to a data file by its index.

        Args:
            index (int or str): Data file index or description.
            filepath (str): Path to the data file.
        """
        if isinstance(index, int) and index in self.data_file_dict:
            if index in self.data_files:
                self.data_files[index] = filepath
        elif isinstance(index, str) and index in self.data_file_dict.values():
            for key, val in self.data_files.items():
                if val == index:
                    self.data_files[key] = filepath
        else:
            raise ValueError(f"index must be an integer or string") 
        
    def get_data_file(self, index):
        """
        Get the path to a data file by its index or description.

        Args:
            index (int or str): Data file index or description.

        Returns:
            str: Path to the data file.
        """
        if isinstance(index, int) and index in self.data_file_dict:
            return self.data_files[index]
        elif isinstance(index, str) and index in self.data_file_dict.values():
            for key, val in self.data_files.items():
                if val == index:
                    return key
        else:
            raise ValueError("index must be an integer or string")
        
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