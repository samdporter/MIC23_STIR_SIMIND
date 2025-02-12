### This file contains functions to convert SIMIND output header files to STIR header files.
### There is some horrible stuff here that I'm not very proud of. I'm sorry.

### Author: Sam Porter

import os
import sys
from pathlib import Path

from sirf.STIR import AcquisitionData


class Converter:
    """
    Class to convert SIMIND header files to STIR header files
    """

    @staticmethod
    def convert_line(line, dir_switch):
        """
        Converts a single line from SIMIND to STIR format
        """

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
        elif "orbit" and "noncircular" in line: 
            res =  f"orbit := non-circular"
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
            # convert to STIR format. If angle +180 > 360, then -360
            angle = dir_switch*float(line.split()[3]) + 180
            if angle >= 360:
                angle -= 360
            res =  f"start angle := {angle}"
        elif "!name of data file" in line:
            file = Path(line.split()[5])
            # only filename plus extension
            res = f"!name of data file := {file.stem+file.suffix}"
        else:
            res = line
        
        return res, dir_switch

    @staticmethod
    def convert(filename, return_object=True):
        """
        Converts a SIMIND header file to a STIR header file
        """

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
        Read two files and align values in file_to_adjust based on reference_file.
        
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


