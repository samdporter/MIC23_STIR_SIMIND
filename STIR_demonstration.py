#!/usr/bin/env python
# coding: utf-8

# In[24]:

from sirf.STIR import (ImageData, AcquisitionData,
                       SPECTUBMatrix, AcquisitionModelUsingMatrix,
                       MessageRedirector,)
from src.simind import *
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd

import argparse

msg = MessageRedirector()
# AcquisitionData.set_storage_scheme('memory')

# In[25]:

parser = argparse.ArgumentParser(description='Run a simulation using SIMIND and STIR')

parser.add_argument('--total_activity', type=float, default=258.423, help='Total activity in MBq')
parser.add_argument('--time_per_projection', type=int, default=43, help='Time per projection in seconds')
parser.add_argument('--photon_multiplier', type=float, default=0.001, help='Number of photons simulated is calculated based on source map. This number multiplies the calculated number of photons')
parser.add_argument('--photopeak_energy', type=float, default=208, help='Photopeak energy in keV')
parser.add_argument('--window_lower', type=float, default=187.56, help='Lower window in keV')
parser.add_argument('--window_upper', type=float, default=229.24, help='Upper window in keV')
parser.add_argument('--source_type', type=str, default='lu177', help='Source type')
parser.add_argument('--collimator', type=str, default='G8-MEGP', help='Collimator')
parser.add_argument('--kev_per_channel', type=float, default=10., help='keV per channel')
parser.add_argument('--max_energy', type=float, default=498.3, help='Max energy in keV')
parser.add_argument('--mu_map_path', type=str, default='data/Lu177/registered_CTAC.hv', help='Path to mu map')
parser.add_argument('--image_path', type=str, default='data/Lu177/osem_reconstruction_postfilter_555.hv', help='Path to image')
parser.add_argument('--measured_data_path', type=str, default='data/Lu177/SPECTCT_NEMA_128_EM001_DS_en_1_Lu177_EM.hdr', help='Path to measured data')
parser.add_argument('--measured_additive', type=str, default='/home/sam/working/STIR_users_MIC2023/data/Lu177/STIR_TEW.hs', help='Path to measured additive')
parser.add_argument('--output_dir', type=str, default='simind_output', help='Output directory')
parser.add_argument('--output_prefix', type=str, default='output', help='Output prefix')
parser.add_argument('--input_smc_file_path', type=str, default='input/input.smc', help='Path to input smc file')
parser.add_argument('--scoring_routine', type=int, default=1, help='Scoring routine')
parser.add_argument('--collimator_routine', type=int, default=1, help='Collimator routine')
parser.add_argument('--photon_direction', type=int, default=3, help='Photon direction')
parser.add_argument('--crystal_thickness', type=float, default=7.25, help='Crystal thickness in mm')
parser.add_argument('--crystal_half_length_radius', type=float, default=393.6/2, help='Crystal half length radius in mm')
parser.add_argument('--crystal_half_width', type=float, default=511.7/2, help='Crystal half width in mm')
parser.add_argument('--flag_11', type=bool, default=True, help='Flag 11 - use collimator')

args = parser.parse_args()

num_energy_spectra_channels = args.max_energy//args.kev_per_channel

def get_acquisition_model(measured_data, additive_data, image, mu_map_stir):
    acq_matrix = SPECTUBMatrix()
    acq_matrix.set_attenuation_image(mu_map_stir)
    acq_matrix.set_keep_all_views_in_cache(True)
    acq_matrix.set_resolution_model(1.81534, 0.02148, False)
    try:
        acq_model.set_additive(additive_data)
    except Exception as e:
        print(e)
        print("Could not set additive data")
    acq_model = AcquisitionModelUsingMatrix(acq_matrix)
    acq_model.set_up(measured_data, image)
    return acq_model

# In[26]:

def main(args):


    image = ImageData(args.image_path)
    mu_map = ImageData(args.mu_map_path)
    measured_data = AcquisitionData(args.measured_data_path)
    measured_additive = AcquisitionData(args.measured_additive)

    # In[27]:

    ## Unfortunately this is necessary due to a bug in STIR
    # Only for the STIR reconstruction
    mu_map_stir = mu_map.clone()
    mu_map_stir.fill(np.flip(mu_map.as_array(), axis=2))

    # In[29]:


    os.chdir("/home/sam/working/STIR_users_MIC2023")
    # set up the simulator
    simulator = SimindSimulator(template_smc_file_path=args.input_smc_file_path,
                                output_dir=args.output_dir, output_prefix= args.output_prefix,
                                source=image, mu_map=mu_map, template_sinogram=measured_data)


    # Set the parameters for the simulation
    simulator.add_comment("Demonstration of SIMIND simulation")
    simulator.set_windows(args.window_lower, args.window_upper, 0)
    simulator.add_index("photon_energy", args.photopeak_energy)
    simulator.add_index("scoring_routine", args.scoring_routine)
    simulator.add_index("collimator_routine", args.collimator_routine)
    simulator.add_index("photon_direction", args.photon_direction)
    simulator.add_index("source_activity", args.total_activity*args.time_per_projection)
    # crystal dimensions
    simulator.add_index("crystal_thickness", args.crystal_thickness/10) # cm
    simulator.add_index("crystal_half_length_radius", args.crystal_half_length_radius/10)
    simulator.add_index("crystal_half_width", args.crystal_half_width/10)
    simulator.config.set_flag(11, args.flag_11)
    simulator.add_index("step_size_photon_path_simulation", min(*image.voxel_sizes())/10) # cm
    # resolutin
    simulator.add_index("energy_resolution", 9.5) # percent
    simulator.add_index("intrinsic_resolution", 0.31) # cm


    # Set the runtime switches
    simulator.add_runtime_switch("CC", args.collimator) # which collimator
    simulator.add_runtime_switch("NN", args.photon_multiplier) # multiplier for number of photons simulated
    simulator.add_runtime_switch("FI", args.source_type) # source type

    simulator.run_simulation()

    simind_total = simulator.get_output_total()
    simind_scatter = simulator.get_output_scatter()
    simind_true = simind_total - simind_scatter

    base_output_filename = f"NN{args.photon_multiplier}_CC{args.collimator}_FI{args.source_type}_"
    simind_total.write(os.path.join(args.output_dir, "simind_total_" + base_output_filename))
    simind_scatter.write(os.path.join(args.output_dir, "simind_scatter_" + base_output_filename))
    simind_true.write(os.path.join(args.output_dir, "simind_true_" + base_output_filename))


    acq_model = get_acquisition_model(measured_data, measured_additive, image, mu_map_stir)
    stir_forward_projection = acq_model.forward(image)

    print(f"simind total counts: {simind_total.sum()}")
    print(f"simind true counts: {simind_true.sum()}")
    print(f"simind scatter counts: {simind_scatter.sum()}")
    print("\n")
    print("\n")
    print(f"measured total counts: {measured_data.sum()}")
    print(f"stir true counts: {stir_forward_projection.sum()}")
    print(f"measured additive counts: {measured_additive.sum()}")

    # In[45]:

    # save counts to csv
    counts = pd.DataFrame({
        "simind_total": [simind_total.sum()],
        "simind_true": [simind_true.sum()],
        "simind_scatter": [simind_scatter.sum()],
        "measured": [measured_data.sum()],
        "stir_forward": [stir_forward_projection.sum()],
        "measured_additive": [measured_additive.sum()]
    })

    counts.to_csv(os.path.join(args.output_dir, base_output_filename + ".csv"))

    data_list = [
        ((simind_total), "simind total"),
        ((simind_true), "simind true"),
        ((simind_scatter), "simind scatter"),
        ((measured_data), "measured"),
        ((stir_forward_projection), "stir forward"),
        ((measured_additive), "measured additive")
    ]

    data_list = [(data.as_array(), title) for data, title in data_list]

    axial_slice = 55

    vmax = max([data[0][axial_slice].max() for data, _ in data_list])


    # In[48]:

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

    for i, (data, title) in enumerate(data_list):
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


    # In[ ]:

if __name__ == '__main__':
    
    try:
        main(args)
    except Exception as e:
        print(e)
        raise e