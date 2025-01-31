#%% 
# This cell imports necessary libraries and sets up the environment for SIMIND and STIR simulation.

from sirf.STIR import (ImageData, AcquisitionData,
                       SPECTUBMatrix, AcquisitionModelUsingMatrix,
                       MessageRedirector, RelativeDifferencePrior,)
from src.simulator import *
from src.simind_projector import *
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np
import os

import argparse

msg = MessageRedirector()

#%% 
def get_args():
    parser = argparse.ArgumentParser(description='Simulation Parameters')

    parser.add_argument('--initial_activity', type=float, default=258.423, help='Initial activity in MBq')
    parser.add_argument('--time_per_projection', type=int, default=43, help='Time per projection in seconds')
    parser.add_argument('--photon_multiplier', type=int, default=10, help='Number of photons simulated')
    parser.add_argument('--photopeak_energy', type=float, default=208, help='Photopeak energy in keV')
    parser.add_argument('--window_lower', type=float, default=187.56, help='Lower window energy in keV')
    parser.add_argument('--window_upper', type=float, default=229.24, help='Upper window energy in keV')
    parser.add_argument('--source_type', type=str, default='lu177', help='Source type')
    parser.add_argument('--collimator', type=str, default='G8-MEGP', help='Collimator type')
    parser.add_argument('--kev_per_channel', type=float, default=5.0, help='Energy per channel in keV')
    parser.add_argument('--max_energy', type=float, default=498.3, help='Maximum energy in keV')
    parser.add_argument('--mu_map_path', type=str, default='data/Lu177/registered_CTAC.hv', help='Path to the mu-map file')
    parser.add_argument('--image_path', type=str, default='data/Lu177/osem_image.hv', help='Path to the image file')
    parser.add_argument('--measured_data_path', type=str, default='data/Lu177/SPECTCT_NEMA_128_EM001_DS_en_1_Lu177_EM.hdr', help='Path to the measured data file')
    parser.add_argument('--measured_additive_path', type=str, default='data/Lu177/tew_scatter.hs', help='Path to the measured additive file')
    parser.add_argument('--output_dir', type=str, default='simind_output', help='Output directory')
    parser.add_argument('--output_prefix', type=str, default='output', help='Output file prefix')
    parser.add_argument('--input_smc_file_path', type=str, default='input/input.smc', help='Path to the input SMC file')
    parser.add_argument('--scoring_routine', type=int, default=1, help='Scoring routine')
    parser.add_argument('--collimator_routine', type=int, default=0, help='Collimator routine (1 for septal penetration)')
    parser.add_argument('--photon_direction', type=int, default=2, help='Photon acceptance angle')
    parser.add_argument('--crystal_thickness', type=float, default=7.25, help='Crystal thickness in mm')
    parser.add_argument('--crystal_half_length_radius', type=float, default=393.6 / 2, help='Crystal half-length radius in mm')
    parser.add_argument('--crystal_half_width', type=float, default=511.7 / 2, help='Crystal half-width in mm')
    parser.add_argument('--flag_11', type=bool, default=True, help='Flag 11 - use collimator')
    parser.add_argument('--num_energy_spectra_channels', type=int, default=int(498.3 // 5.0), help='Number of energy spectra channels')
    parser.add_argument('--num_subsets', type=int, default=16, help='Number of subsets')
    parser.add_argument('--num_epochs', type=int, default=15, help='Number of epochs')
    parser.add_argument('--half_life', type=float, default=6.647*24, help='Half life of the isotope in hours')

    return parser.parse_args()

#%% 
# Helper function to get the acquisition model based on inputs.

def get_acquisition_model(measured_data, additive_data, image, mu_map_stir):
    acq_matrix = SPECTUBMatrix()
    acq_matrix.set_attenuation_image(mu_map_stir)
    acq_matrix.set_keep_all_views_in_cache(True)
    #acq_matrix.set_resolution_model(1.81534, 0.02148, False)
    acq_matrix.set_resolution_model(0.9, 0.02, False)
    acq_model = AcquisitionModelUsingMatrix(acq_matrix)
    try:
        acq_model.set_additive_term(additive_data)
    except Exception as e:
        print(e)
        print("Could not set additive data")
    acq_model.set_up(measured_data, image)
    return acq_model

def get_average_activity(measured_data, time_per_projection, initial_activity, half_life):
    
    decay_constant = np.log(2) / half_life
    total_activity = initial_activity / decay_constant * (1 - np.exp(-decay_constant * time_per_projection * measured_data.dimensions()[2]))
    average_activity = total_activity / time_per_projection / measured_data.dimensions()[2]
    
    return average_activity

def save_estimates_tofile(projector, output_dir, suffix):
    """
    Save the estimates in the output directory.
    """
    projector.additive_correction.write(f"{output_dir}/additive_correction_{suffix}.hs")
    projector.additive_estimate.write(f"{output_dir}/additive_estimate_{suffix}.hs")
    
def save_images_tofile(projector, output_dir, suffix):
    """
    Save the images in the output directory.
    """
    projector.simind_simulator.source.write(f"{output_dir}/image_{suffix}.hv")
    projector.simind_simulator.mu_map.write(f"{output_dir}/mu_map_{suffix}.hv")

def calculate_sensitivity(projector, subset, inv_sensitivity_ims, epoch):
    """
    Calculate and store sensitivity images if not already computed.
    """
    if epoch == 0:
        print(f"Calculating sensitivity for subset {subset}")
        inv_sensitivity = projector.backward(projector.acquisition_data.get_uniform_copy(1), subset)
        inv_sensitivity.fill(np.reciprocal(inv_sensitivity.as_array(), where=inv_sensitivity.as_array() != 0))
        inv_sensitivity = inv_sensitivity.maximum(0)
        inv_sensitivity_ims[subset] = inv_sensitivity
    return inv_sensitivity_ims[subset]

def update_estimate(current_estimate, projector, inv_sensitivity, subset, eps):
    """
    Update the current estimate based on the acquisition data and sensitivity.
    """
    acq_ratio = projector.acquisition_data / (projector.forward(current_estimate, subset) + eps)
    backproj_of_acq_ratio = projector.backward(acq_ratio, subset)
    current_estimate = (current_estimate+eps) * backproj_of_acq_ratio * inv_sensitivity
    return current_estimate.maximum(0)

def OSEM(init_image, projector, num_subsets, num_epochs, eps=1e-6, update_interval = 5, 
        residual_correction=False, with_scatter=False, save_estimates=False, save_images=False, output_dir="output"):
    """
    Run OSEM with optional SIMIND corrections.
    """
    current_estimate = init_image.clone()
    inv_sensitivity_ims = [0]*num_subsets

    for j in range(num_epochs):
        for i in range(num_subsets):
            inv_sensitivity = calculate_sensitivity(projector, i, inv_sensitivity_ims, j)
            current_estimate = update_estimate(current_estimate, projector, inv_sensitivity, i, eps)

        if (residual_correction or with_scatter) and j != num_epochs - 1 and j % update_interval == 0 and j != 0:
            if residual_correction:
                projector.residual_correction = True
            if with_scatter:
                projector.update_scatter = True

            projector.simulate_forward_projection(current_estimate.copy())
            if save_estimates:
                save_estimates_tofile(projector, output_dir, f"epoch_{j}")
            if save_images:
                save_images_tofile(projector, output_dir, f"epoch_{j}")

    return current_estimate


#%% 
# Loading the images and data for the simulation.
def main(args):

    image = ImageData(args.image_path)
    mu_map = ImageData(args.mu_map_path)
    measured_data = AcquisitionData(args.measured_data_path)
    measured_additive = AcquisitionData(args.measured_additive_path)

    average_activity = get_average_activity(measured_data, args.time_per_projection, args.initial_activity, args.half_life * 3600)
    print(f"Average activity: {average_activity} MBq")

    #%% 
    # Flipping mu_map along the z-axis due to a bug in STIR.

    mu_map_stir = mu_map.clone()
    mu_map_stir.fill(np.flip(mu_map.as_array(), axis=2))

    #%% 
    # Setting up the SIMIND simulator.

    os.chdir("/home/sam/working/STIR_users_MIC2023")
    simulator = SimindSimulator(template_smc_file_path=args.input_smc_file_path,
                                output_dir=args.output_dir, output_prefix=args.output_prefix,
                                source=image, mu_map=mu_map, template_sinogram=measured_data)

    # Setting up the simulator parameters.

    simulator.add_comment("Demonstration of SIMIND simulation")
    simulator.set_windows(args.window_lower, args.window_upper, 0)
    simulator.add_index("photon_energy", args.photopeak_energy)
    simulator.add_index("scoring_routine", args.scoring_routine)
    simulator.add_index("collimator_routine", args.collimator_routine)
    simulator.add_index("photon_direction", args.photon_direction)
    simulator.add_index("source_activity", args.initial_activity * args.time_per_projection)
    simulator.add_index("crystal_thickness", args.crystal_thickness / 10)  # cm
    simulator.add_index("crystal_half_length_radius", args.crystal_half_length_radius / 10)  # cm
    simulator.add_index("crystal_half_width", args.crystal_half_width / 10)  # cm
    simulator.config.set_flag(11, args.flag_11)
    simulator.add_index("step_size_photon_path_simulation", min(*image.voxel_sizes()) / 10)  # cm
    simulator.add_index("energy_resolution", 9.5)  # percent
    simulator.add_index("intrinsic_resolution", 0.31)  # cm
    simulator.add_index("cutoff_energy_terminate_photon_history", args.window_lower*0.5)  # keV

    simulator.add_runtime_switch("CC", args.collimator)  # which collimator
    simulator.add_text_variable(1, args.collimator) # changes collimator in the .smc file (probably a duplication of the above)
    simulator.add_runtime_switch("NN", args.photon_multiplier)  # photon multiplier
    simulator.add_runtime_switch("FI", args.source_type)  # source type

    #%%
    projector = SimindProjector()
    projector.simind_simulator = simulator
    projector.image = image
    projector.acquisition_data = measured_data

    #%% 
    # Setting up the acquisition model for forward projection and displaying some results.

    acq_model = get_acquisition_model(measured_data, measured_additive, image, mu_map_stir)
    projector.stir_projector = acq_model
    #%% # 

    # %%
    projector.num_subsets = args.num_subsets
    init_image = image.get_uniform_copy(1)
    fov_filter = TruncateToCylinderProcessor()
    fov_filter.apply(init_image)

    #%%
    osem_image = OSEM(init_image, projector, num_subsets=args.num_subsets, num_epochs=args.num_epochs)
    osem_image.write(f"{args.output_dir}/osem_image.hv")
    osem_image = ImageData(f"{args.output_dir}/osem_image.hv")

    #%%
    osem_simind_correction_image = OSEM(init_image, projector, num_subsets=args.num_subsets, num_epochs=args.num_epochs, 
                                        residual_correction=True, 
                                        save_estimates=True, output_dir=args.output_dir)
    osem_simind_correction_image.write(f"{args.output_dir}/osem_simind_correction_image.hv")

    #%%
    #osem_simind_scatter_image = OSEM(init_image, projector, num_subsets=args.num_subsets, num_epochs=args.num_epochs, 
    #                                 with_scatter=True)
    #osem_simind_scatter_image.write(f"{args.output_dir}/osem_simind_scatter_image.hv")

    #osem_simind_scatter_correction_image = OSEM(init_image, projector, num_subsets=args.num_subsets, num_epochs=args.num_epochs,
    #                                            residual_correction=True, with_scatter=True, )
    #osem_simind_scatter_correction_image.write(f"{args.output_dir}/osem_simind_scatter_correction_image.hv")
# %%
# save the images as .hv files
# %%
# display the images

if __name__ == "__main__":

    args = get_args()

    main(args)

