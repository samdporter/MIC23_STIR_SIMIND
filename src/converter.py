import os
import re
import sys
from pathlib import Path
from sirf.STIR import AcquisitionData
import fileinput
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

class Converter:
    """
    Class to convert SIMIND header files to STIR header files.
    """

    @staticmethod
    def convert_line(line, dir_switch):
        """
        Converts a single line from SIMIND to STIR format.
        """
        patterns = {
            "program", "patient", "institution", "contact", "ID", "exam type",
            "detector head", "number of images/energy window", "time per projection",
            "data description", "total number of images", "acquisition mode"
        }
        
        if any(pattern in line for pattern in patterns):
            return ";" + line, dir_switch

        if "Radius" in line:
            return f"Radius := {float(line.split()[-1])}", dir_switch
        elif "orbit" in line and "noncircular" in line:
            return "orbit := non-circular", dir_switch
        elif "!number format := short float" in line:
            return "!number format := float", dir_switch
        elif "image duration" in line:
            parts = line.split()
            return f"number of time frames := 1\nimage duration (sec) [1] := {parts[4]}", dir_switch
        elif ";energy window lower level" in line:
            return f"energy window lower level[1] := {line.split()[-1]}", dir_switch
        elif ";energy window upper level" in line:
            return f"energy window upper level[1] := {line.split()[-1]}", dir_switch
        elif "CCW" in line:
            return line, -1
        elif "start angle" in line:
            angle = dir_switch * float(line.split()[3]) + 180
            return f"start angle := {angle % 360}", dir_switch
        elif "!name of data file" in line:
            file = Path(line.split()[5])
            return f"!name of data file := {file.stem + file.suffix}", dir_switch
        
        return line, dir_switch

    @staticmethod
    def convert(filename, return_object=True):
        """
        Converts a SIMIND header file to a STIR header file.
        """
        if not filename.endswith(".h00"):
            logging.error("USAGE: script filename.h00")
            sys.exit(1)

        stirfilename = filename.replace(".h00", ".hs")
        dir_switch = 1

        with open(filename, "r") as f_in, open(stirfilename, "w") as f_out:
            for line in f_in:
                write_line, dir_switch = Converter.convert_line(line.strip(), dir_switch)
                f_out.write(write_line + "\n")

        logging.info(f"Output written to {stirfilename}")
        return AcquisitionData(stirfilename) if return_object else None

    @staticmethod
    def adjust_values(reference_file, file_to_adjust, threshold=None, output_adjusted_file=None):
        """
        Adjust values in file_to_adjust based on reference_file, either replacing exactly or within a threshold.
        """
        if output_adjusted_file is None:
            output_adjusted_file = file_to_adjust[:-4] + "_adjusted.h00"
            
        if isinstance(reference_file, AcquisitionData):
            reference_file.write("tmp_ref.hs")
            reference_file = "tmp_ref.hs"

        with open(reference_file, "r") as ref_file, open(file_to_adjust, "r") as to_adjust_file:
            reference_lines = {line.split(":=")[0].strip(): line.split(":=")[1].strip() for line in ref_file if ":=" in line}
            adjust_lines = to_adjust_file.readlines()
        
        for i, line in enumerate(adjust_lines):
            if ":=" in line:
                key, value = map(str.strip, line.split(":="))
                if key in reference_lines:
                    try:
                        ref_value, adj_value = float(reference_lines[key]), float(value)
                        if threshold is None or abs(ref_value - adj_value) <= threshold:
                            adjust_lines[i] = f"{key} := {ref_value}\n"
                    except ValueError:
                        pass
        
        with open(output_adjusted_file, "w") as out:
            out.writelines(adjust_lines)
            
        if isinstance(reference_file, str) and "tmp_ref.hs" in reference_file:
            os.remove(reference_file)
        
        return AcquisitionData(output_adjusted_file)

    @staticmethod
    def replace_sinogram_values(reference_sinogram, sinogram_to_adjust):
        """
        Replace values in a sinogram file based on a reference sinogram.
        """
        ref_filename, adj_filename = "tmp_ref.hs", "tmp_adjust.hs"
        reference_sinogram.write(ref_filename)
        sinogram_to_adjust.write(adj_filename)
        result = Converter.adjust_values(ref_filename, adj_filename)
        for tmp in [ref_filename, adj_filename, "tmp_ref.s", "tmp_adjust.s"]:
            os.remove(tmp)
        return result

    @staticmethod
    def convert_sinogram_parameter(sinogram, parameter, value):
        """
        Modify a sinogram parameter in-place.
        """
        filename = "tmp.hs"
        sinogram.write(filename)
        pattern = re.compile(rf"^!?{re.escape(parameter)}\s*:=", re.IGNORECASE)
        
        with fileinput.input(filename, inplace=True) as file:
            for line in file:
                if pattern.match(line.strip()):
                    print(f"{parameter} := {value}")
                else:
                    print(line, end="")
        
        sinogram = AcquisitionData(filename)
        os.remove(filename)
        return sinogram