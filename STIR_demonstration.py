#!/usr/bin/env python
# coding: utf-8

from sirf.STIR import (ImageData, AcquisitionData,
                       SPECTUBMatrix, AcquisitionModelUsingMatrix,
                       MessageRedirector,)
from src.simulator import SimindSimulator
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
import numpy as np
import os

import argparse

msg = MessageRedirector()
# AcquisitionData.set_storage_scheme('memory')

parser = argparse.ArgumentParser(description='Run a simulation using SIMIND and STIR')

parser.add_argument('--total_activity', type=float, default=258.423, help='Total activity in MBq')
parser.add_argument('--time_per_projection', type=int, default=43, help='Time per projection in seconds')
parser.add_argument('--photon_multiplier', type=float, default=10., help='Number of photons simulated is calculated based on source map. This number multiplies the calculated number of photons')
parser.add_argument('--photopeak_energy', type=float, default=208, help='Photopeak energy in keV')
parser.add_argument('--window_lower', type=float, default=187.56, help='Lower window in keV')
parser.add_argument('--window_upper', type=float, default=229.24, help='Upper window in keV')
parser.add_argument('--source_type', type=str, default='lu177', help='Source type')
parser.add_argument('--collimator', type=str, default='G8-MEGP', help='Collimator')
parser.add_argument('--kev_per_channel', type=float, default=10., help='keV per channel')
parser.add_argument('--max_energy', type=float, default=498.3, help='Max energy in keV')
parser.add_argument('--mu_map_path', type=str, default='data/Lu177/registered_CTAC.hv', help='Path to mu map')
parser.add_argument('--image_path', type=str, default='data/Lu177/osem_image.hv', help='Path to image')
parser.add_argument('--measured_data_path', type=str, default='data/Lu177/SPECTCT_NEMA_128_EM001_DS_en_1_Lu177_EM.hdr', help='Path to measured data')
parser.add_argument('--measured_additive_path', type=str, default='data/Lu177/tew_scatter.hs', help='Path to measured additive')
parser.add_argument('--output_dir', type=str, default='simind_output', help='Output directory')
parser.add_argument('--output_prefix', type=str, default='output', help='Output prefix')
parser.add_argument('--input_smc_file_path', type=str, default='input/input.smc', help='Path to input smc file')
parser.add_argument('--scoring_routine', type=int, default=1, help='Scoring routine')
parser.add_argument('--collimator_routine', type=int, default=0, help='Collimator routine')
parser.add_argument('--photon_direction', type=int, default=3, help='Photon direction')
parser.add_argument('--crystal_thickness', type=float, default=7.25, help='Crystal thickness in mm')
parser.add_argument('--crystal_half_length_radius', type=float, default=393.6/2, help='Crystal half length radius in mm')
parser.add_argument('--crystal_half_width', type=float, default=511.7/2, help='Crystal half width in mm')
parser.add_argument('--flag_11', type=bool, default=True, help='Flag 11 - use collimator')
parser.add_argument('--half_life', type=float, default=6.647*24, help='Half life of the isotope in hours')
parser.add_argument('--axial_slice', type=int, default=65, help='Axial slice to plot')

args = parser.parse_args()

num_energy_spectra_channels = args.max_energy//args.kev_per_channel

def get_acquisition_model(measured_data, additive_data, image, mu_map_stir):
    acq_matrix = SPECTUBMatrix()
    acq_matrix.set_attenuation_image(mu_map_stir)
    acq_matrix.set_keep_all_views_in_cache(True)
    acq_matrix.set_resolution_model(0.9323, 0.03, False) 
    acq_model = AcquisitionModelUsingMatrix(acq_matrix)
    try:
        acq_model.set_additive_term(additive_data)
    except Exception as e:
        print(e)
        print("Could not set additive data")
    acq_model.set_up(measured_data, image)
    return acq_model

def lower_threshold_image(image, threshold):
    """ saves lots of time in simulation"""
    image_array = image.as_array()
    image_array[image_array < threshold] = 0
    image.fill(image_array)
    return image

def main(args):
    image = ImageData(args.image_path)
    # threshold to 1% of max value
    image = lower_threshold_image(image, 0.01*image.max())
    mu_map = ImageData(args.mu_map_path)
    measured_data = AcquisitionData(args.measured_data_path)
    
    measured_additive = None
    if args.measured_additive_path and args.measured_additive_path.endswith(".hs"):
        measured_additive = AcquisitionData(args.measured_additive_path)

    os.chdir("/home/sam/working/STIR_users_MIC2023")

    simulator = SimindSimulator(
        template_smc_file_path=args.input_smc_file_path,
        output_dir=args.output_dir, 
        output_prefix=args.output_prefix,
        source=image, 
        mu_map=mu_map, 
        template_sinogram=measured_data
    )

    simulator.add_comment("Demonstration of SIMIND simulation")
    simulator.set_windows(args.window_lower, args.window_upper, 0)
    simulator.add_index("photon_energy", args.photopeak_energy)
    simulator.add_index("scoring_routine", args.scoring_routine)
    simulator.add_index("collimator_routine", args.collimator_routine)
    simulator.add_index("photon_direction", args.photon_direction)
    simulator.add_index("source_activity", args.total_activity * args.time_per_projection)
    simulator.add_index("crystal_thickness", args.crystal_thickness / 10) 
    simulator.add_index("crystal_half_length_radius", args.crystal_half_length_radius / 10)
    simulator.add_index("crystal_half_width", args.crystal_half_width / 10)
    simulator.config.set_flag(11, args.flag_11)
    simulator.add_index("step_size_photon_path_simulation", min(*image.voxel_sizes()) / 10) 
    simulator.add_index("energy_resolution", 9.5)
    simulator.add_index("intrinsic_resolution", 0.31)
    simulator.add_index("cutoff_energy_terminate_photon_history", args.window_lower * 0.5)
    
    simulator.add_runtime_switch("CC", args.collimator)
    simulator.add_runtime_switch("NN", args.photon_multiplier)
    simulator.add_runtime_switch("FI", args.source_type)
    
    simulator.run_simulation()

    simind_total = simulator.get_output_total()
    simind_scatter = simulator.get_output_scatter()
    simind_true = simind_total - simind_scatter

    base_output_filename = f"NN{args.photon_multiplier}_CC{args.collimator}_FI{args.source_type}_"

    counts = {
        "simind_total": simind_total.sum(),
        "simind_true": simind_true.sum(),
        "simind_scatter": simind_scatter.sum(),
    }
    
    if measured_additive is not None:
        counts["measured_additive"] = measured_additive.sum()
    if args.measured_additive_path and args.measured_additive_path.endswith(".hs"):
        counts["measured"] = measured_data.sum()
    
    pd.DataFrame([counts]).to_csv(os.path.join(args.output_dir, base_output_filename + ".csv"))

    data_list = [
        (simind_total, "simind total"),
        (measured_data, "measured"),
        (simind_true, "simind true"),
        (measured_additive.maximum(0) if measured_additive is not None else None, "measured scatter"),
        (simind_scatter, "simind scatter"),
        (measured_additive if measured_additive is not None else None, "TEW scatter"),
    ]
    
    # Filter out None values
    data_list = [(data.as_array(), title) for data, title in data_list if data is not None]

    axial_slice = args.axial_slice
    vmax = max([data[0][axial_slice].max() for data, _ in data_list])

    # Define consistent font size and colormap
    font_size = 14
    colormap = 'viridis'

    # Create a figure and a GridSpec with 3 rows
    fig = plt.figure(figsize=(len(data_list)*4,7*2,))
    gs = GridSpec(3, len(data_list), height_ratios=[2, 0.15, 3])  # Adjusted GridSpec for clarity

    # Create image subplots in the first row
    ax_images = [fig.add_subplot(gs[0, i]) for i in range(len(data_list))]

    for i, (data, title) in enumerate(data_list):
        im = ax_images[i].imshow(data[0, axial_slice], vmin=0, vmax=vmax, cmap=colormap)
        ax_images[i].set_title(f"{title}: {np.trunc(data.sum())} ", fontsize=font_size)
        ax_images[i].axis('off')

    # Place a colorbar in a new row, just for the colorbar
    cbar_ax = fig.add_subplot(gs[1, :])  # Spanning across the bottom of the image plots
    fig.colorbar(im, cax=cbar_ax, orientation='horizontal', pad=0.02)  # Reduced padding
    cbar_ax.set_xlabel('Counts', fontsize=font_size)
    cbar_ax.xaxis.set_label_position('top')

    # Set consistent font size and line width for the line plot
    line_width = 2

    # Create line plot in the third row
    ax_line = fig.add_subplot(gs[2, :])  # Spanning across both columns in the third row

    # Plotting
    # make colours eqwually spread out with length of data_list
    colours = plt.cm.viridis(np.linspace(0, 1, len(data_list)))
    for i, (data, title) in enumerate(data_list):
        #ax_line.plot(data[40], label=title, linewidth=line_width, color=colours[i], linestyle='--')
        ax_line.plot(data[0, axial_slice][60], linewidth=line_width, color=colours[i], linestyle='-', label=title)
        #ax_line.plot(data[90], linewidth=line_width, color=colours[i], linestyle=':')

    # Enhance the appearance of the line plot
    ax_line.set_xlabel('Projection angle', fontsize=font_size)
    ax_line.set_ylabel('Intensity', fontsize=font_size)
    ax_line.set_title(f'Profile Through Sinogram', fontsize=font_size + 2)
    ax_line.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax_line.legend(loc='upper left', fontsize=font_size)
    # et minimum and maximum values for x-axis
    ax_line.set_xlim(0, 128)

    # Adjust spacing and layout
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "comparison_" + base_output_filename + ".png"))
    plt.close()

    # In[ ]:

    coronal_slice = 55
    vmax = max([data[0][:, coronal_slice].max() for data, _ in data_list])

    # Create a figure and a GridSpec with 3 rows
    fig = plt.figure(figsize=(len(data_list)*4,7*2,))
    gs = GridSpec(3, len(data_list), height_ratios=[2, 0.15, 3])  # Adjusted GridSpec for clarity

    # Create image subplots in the first row
    ax_images = [fig.add_subplot(gs[0, i]) for i in range(len(data_list))]

    filtered_data = [(data, title) for data, title in data_list if data is not None]
    for i, (data, title) in enumerate(filtered_data):
        im = ax_images[i].imshow(data[0, :, coronal_slice], vmin=0, vmax=vmax, cmap=colormap)
        ax_images[i].set_title(f"{title}: {np.trunc(data.sum())} ", fontsize=font_size)
        ax_images[i].axis('off')

    # Place a colorbar in a new row, just for the colorbar
    cbar_ax = fig.add_subplot(gs[1, :])  # Spanning across the bottom of the image plots
    fig.colorbar(im, cax=cbar_ax, orientation='horizontal', pad=0.02)  # Reduced padding
    cbar_ax.set_xlabel('Counts', fontsize=font_size)
    cbar_ax.xaxis.set_label_position('top')

    # Set consistent font size and line width for the line plot
    line_width = 2

    # Create line plot in the third row
    ax_line = fig.add_subplot(gs[2, :])  # Spanning across both columns in the third row

    # Plotting
    # make colours eqwually spread out with length of data_list
    colours = plt.cm.viridis(np.linspace(0, 1, len(data_list)))
    for i, (data, title) in enumerate(data_list):
        #ax_line.plot(data[40], label=title, linewidth=line_width, color=colours[i], linestyle='--')
        ax_line.plot(data[0, 60, coronal_slice], linewidth=line_width, color=colours[i], linestyle='-', label=title)
        #ax_line.plot(data[90], linewidth=line_width, color=colours[i], linestyle=':')

    # Enhance the appearance of the line plot
    ax_line.set_xlabel('Projection angle', fontsize=font_size)
    ax_line.set_ylabel('Intensity', fontsize=font_size)
    ax_line.set_title(f'Profile Through Sinogram', fontsize=font_size + 2)
    ax_line.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax_line.legend(loc='upper left', fontsize=font_size)
    # et minimum and maximum values for x-axis
    ax_line.set_xlim(0, 128)

    # Adjust spacing and layout
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "comparison_coronal_" + base_output_filename + ".png"))
    plt.close()

if __name__ == '__main__':
    
    try:
        main(args)
    except Exception as e:
        print(e)
        raise e