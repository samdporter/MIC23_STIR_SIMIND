import pytest
from sirf.STIR import AcquisitionData

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from functions import extract_attributes_from_stir, STIR_ATTRIBUTE_MAPPING

def compare_stir_attributes(attributes_path, attributes_sinogram):
    """
    Compare STIR attributes extracted from a header file and a sinogram object.
    
    Args:
        attributes_path (dict): Extracted attributes from the header file.
        attributes_sinogram (dict): Extracted attributes from the sinogram object.
    
    Raises:
        AssertionError if any mapped attribute does not match.
    """
    mapped_keys = set(STIR_ATTRIBUTE_MAPPING.keys()) | set(STIR_ATTRIBUTE_MAPPING.values())
    common_keys = mapped_keys & set(attributes_path.keys()) & set(attributes_sinogram.keys())

    for key in common_keys:
        mapped_key = STIR_ATTRIBUTE_MAPPING.get(key, key)  # Use mapping if available
        value_path = attributes_path[key]
        value_sinogram = attributes_sinogram[mapped_key]

        if isinstance(value_path, float) or isinstance(value_sinogram, float):
            assert value_path == pytest.approx(value_sinogram, rel=1e-3), f"Mismatch in {mapped_key}: {value_path} != {value_sinogram}"
        else:
            assert value_path == value_sinogram, f"Mismatch in {mapped_key}: {value_path} != {value_sinogram}"

def test_extract_attributes_from_stir():
    """
    Test that attributes extracted from a STIR header file match those extracted from a sinogram.
    """
    path = "/home/sam/working/STIR_users_MIC2023/data/Lu177/SPECTCT_NEMA_128_EM001_DS_en_1_Lu177_EM.hdr"
    
    attributes_path = extract_attributes_from_stir(path)
    attributes_sinogram = extract_attributes_from_stir(AcquisitionData(path))

    compare_stir_attributes(attributes_path, attributes_sinogram)