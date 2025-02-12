### This file contains functions to convert Hounsfield Units (HU) to attenuation coefficients and densities.
### The functions are based on the bilinear model for attenuation coefficients and densities.
### I wouldn't necessarily trust all this implicitly, but it's a good starting point.
### STIR has it's own way of doing this so TODO might be to wrap that. Or a comparison of the two.

### Author: Sam Porter

import numpy as np
import os

def interpolate_attenuation_coefficient(filename, energy):
    # Skip the header lines; use the 'skiprows' parameter
    energies, coeffs = np.loadtxt(filename, unpack=True, skiprows=12)
    return np.interp(energy, energies, coeffs)

def get_attenuation_coefficient(material, energy, file_path = ''):

    density_water = 1.0  # g/cm^3
    density_bone = 1.85  # g/cm^3 for cortical bone, this is a common average

    if material == 'water':
        filename = os.path.join(file_path, 'h2o.atn')
    elif material == 'bone':
        filename = os.path.join(file_path, 'bone.atn')
    else:
        raise ValueError("Unknown material. Accepted values are 'water' or 'bone'.")
    
    mass_attn_coeffs =  interpolate_attenuation_coefficient(filename, energy)
    if material == 'water':
        return mass_attn_coeffs * density_water
    elif material == 'bone':
        return mass_attn_coeffs * density_bone

def hu_to_attenuation(image_array, photon_energy, file_path = ''):
    # Constants
    HU_water = 0  # HU for water
    HU_bone = 1000  # Typical HU for bone (can vary)

    #convert photon_energy to MeV from keV
    photon_energy = photon_energy / 1000
    
    # Get attenuation coefficients
    mu_water = get_attenuation_coefficient('water', photon_energy, file_path)
    mu_bone = get_attenuation_coefficient('bone', photon_energy, file_path)

    print(mu_water, mu_bone)

    slope_soft = mu_water / 1000  # from -1000 HU (air) to 0 HU (water)
    slope_bone = (mu_bone - mu_water) / 1000  # from 0 HU (water) to 1000 HU (bone)

    # Compute attenuation map using the bilinear model
    attenuation_map = np.where(image_array <= HU_water, mu_water + slope_soft * image_array, mu_water + slope_bone * (image_array - HU_water))
    return attenuation_map

def hu_to_density(image_array):
    # Constants (Typical values; may need to be refined)
    density_air = 0.001225  # g/cm^3
    density_water = 1.0  # g/cm^3
    density_bone = 1.85  # g/cm^3

    slope_soft = (density_water - density_air) / 1000  # from -1000 HU (air) to 0 HU (water)
    slope_bone = (density_bone - density_water) / 1000  # from 0 HU (water) to 1000 HU (bone)

    density_map = np.where(image_array <= 0, density_water + slope_soft * image_array, density_water + slope_bone * (image_array - 0))
    return density_map

def attenuation_to_density(attenuation_array, photon_energy, file_path = ''):
    # Convert photon_energy to MeV from keV
    photon_energy = photon_energy / 1000

    # Obtain attenuation coefficients for water and bone
    mu_water = get_attenuation_coefficient('water', photon_energy, file_path)
    mu_bone = get_attenuation_coefficient('bone', photon_energy, file_path)

    # Densities
    density_air = 0.001225  # g/cm^3
    density_water = 1.0  # g/cm^3
    density_bone = 1.85  # g/cm^3 for cortical bone, this is a common average

    # Slopes for bilinear model
    slope_soft = (density_water - density_air) / mu_water  # from air to water
    slope_bone = (density_bone - density_water) / (mu_bone - mu_water)  # from water to bone

    density_map = np.where(attenuation_array <= mu_water, density_air + slope_soft * attenuation_array, 
                           density_water + slope_bone * (attenuation_array - mu_water))
    return density_map