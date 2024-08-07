#!/usr/bin/env python
# coding: utf-8

# In[46]:

from sirf.STIR import (ImageData, AcquisitionData,
                       SPECTUBMatrix, AcquisitionModelUsingMatrix,
                       MessageRedirector,)
from simind import *
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import argparse

parser = argparse.ArgumentParser(description='Run a simulation of a Y90 SPECT acquisition')
parser.add_argument('--output', type=str, default="data/Y90/", help='Output directory')
parser.add_argument('--source', type=str, default="/home/sam/data/phantom_data/SPECT", help='Source directory')
parser.add_argument('--output_filename', type=str, default="ellipses_megp_cpd", help='Output filename')
parser.add_argument('--image_filename', type=str, default="ellipsoid_image_s.hv", help='Image filename')
parser.add_argument('--umap_filename', type=str, default="umap_zoomed.hv", help='uMap filename')
parser.add_argument('--measured_filename', type=str, default="peak_1_projdata__f1g1d0b0.hs", help='Measured data filename')

# SIMIND parameters
parser.add_argument('--total_activity', type=float, default=187, help='Total activity in MBq')
parser.add_argument('--time_per_projection', type=int, default=20, help='Time per projection in seconds')
parser.add_argument('--photon_multiplier', type=int, default=10, help='Photon multiplier')
parser.add_argument('--photopeak_energy', type=int, default=150, help='Photopeak energy in keV')
parser.add_argument('--window_lower', type=int, default=75, help='Lower window in keV')
parser.add_argument('--window_upper', type=int, default=225, help='Upper window in keV')
parser.add_argument('--source_type', type=str, default="y90", help='Source type')
parser.add_argument('--collimator', type=str, default='ma-megp', help='Collimator')
parser.add_argument('--kev_per_channel', type=int, default=5, help='keV per channel')
parser.add_argument('--num_energy_spectra_channels', type=int, default=400, help='Number of energy spectra channels')

args = parser.parse_args()

# msg = MessageRedirector()
# AcquisitionData.set_storage_scheme('memory')

# In[47]:

# SIMIND parameters
total_activity = args.total_activity  # MBq
time_per_projection = args.time_per_projection  # seconds
photon_multiplier = args.photon_multiplier
photopeak_energy = args.photopeak_energy  # keV
window_lower = args.window_lower  # keV
window_upper = args.window_upper  # keV
source_type = args.source_type
collimator = args.collimator
kev_per_channel = args.kev_per_channel  # keV
num_energy_spectra_channels = args.num_energy_spectra_channels

# In[48]:
image = ImageData(os.path.join(args.source, args.image_filename))
mu_map = ImageData(os.path.join(args.source, args.umap_filename))
measured_data = AcquisitionData(os.path.join(args.source, args.measured_filename))

# In[49]:

mu_map_stir = mu_map.clone()
mu_map_stir.fill(np.flip(mu_map.as_array(), axis=2))

# In[50]:

simulator = SimindSimulator(input_filepath=".", output_filepath=".",)
attributes = {"SourceMap": image,
              "MuMap": mu_map,
              "keVPerChannel": kev_per_channel, 
              'NumberOfEnergySpectraChannels': num_energy_spectra_channels,
              "PhotopeakEnergy": photopeak_energy,
              "SourceActivity": total_activity,
              "PhotonMultiplier": photon_multiplier,
              "ImageDurationPerProjection": time_per_projection,
              "SourceType": source_type.lower(),
              "Collimator": collimator,}
simulator.set_template_sinogram(measured_data)
simulator.set_attributes(attributes)
simulator.set_windows(window_lower, window_upper, 0)

# In[51]:

simulator.run_simulation()

# In[65]:

simind_total = simulator.get_output_total()
simind_scatter = simulator.get_output_scatter()
simind_true = simind_total - simind_scatter

simind_total.write(os.path.join(args.output, f"{args.output_filename}_total.hs"))
simind_scatter.write(os.path.join(args.output, f"{args.output_filename}_scatter.hs"))
simind_true.write(os.path.join(args.output, f"{args.output_filename}_true.hs"))
                  
# In[53]:
acq_matrix = SPECTUBMatrix()
acq_matrix.set_attenuation_image(mu_map_stir)
acq_matrix.set_keep_all_views_in_cache(True)
acq_matrix.set_resolution_model(1.81534, 0.02148, False)
acq_model = AcquisitionModelUsingMatrix(acq_matrix)
acq_model.set_up(measured_data, image)

# In[54]:
stir_forward_projection = acq_model.forward(image)

# In[55]:

scaling_factor_stir = stir_forward_projection.sum()/simind_true.sum()
scaling_factor_measured = measured_data.sum()/simind_total.sum()
scaling_factor = scaling_factor_measured

# In[56]:

print(f"stir scaling factor: {scaling_factor_stir}")
print(f"measured scaling factor: {scaling_factor_measured}")

# In[57]:

print(f"simind unnormalised total counts: {simind_total.sum()}")
print(f"simind unnormalised true counts: {simind_true.sum()}")
print(f"simind unnormalised scatter counts: {simind_scatter.sum()}")
print("\n")
print(f"simind normalised total counts: {simind_total.sum()*scaling_factor}")
print(f"simind normalised true counts: {simind_true.sum()*scaling_factor}")
print(f"simind normalised scatter counts: {simind_scatter.sum()*scaling_factor}")
print("\n")
print(f"measured total counts: {measured_data.sum()}")
print(f"stir total counts: {stir_forward_projection.sum()}")

# In[64]:

# Define consistent font size and colormap
font_size = 14
colormap = 'viridis'
axial_slice = 55

# Set the maximum intensity for color normalization
vmax = 100#max(measured_data.max(), simind_total.max(), stir_forward_projection.max())

data_list = [
    ((simind_total).as_array(), "SIMIND total"),
    (measured_data.as_array(), "measured"),
    ((simind_true).as_array(), "SIMIND true"),
    ((stir_forward_projection).as_array(), "STIR forward"),
    ((simind_scatter).as_array(), "SIMIND scatter"),
]

divisors = [1/scaling_factor_measured, 1, 1/scaling_factor_measured, 1/scaling_factor_stir, 1/scaling_factor_measured]

# Create a figure and a GridSpec with 3 rows
fig = plt.figure(figsize=(len(data_list)*4,7*2,))
gs = GridSpec(3, len(data_list), height_ratios=[2, 0.15, 3])  # Adjusted GridSpec for clarity

# Create image subplots in the first row
ax_images = [fig.add_subplot(gs[0, i]) for i in range(len(data_list))]

for i, (data, title) in enumerate(data_list):
    im = ax_images[i].imshow(data[0, axial_slice]/divisors[i], vmin=0, vmax=vmax, cmap=colormap)
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
    ax_line.plot(data[0, axial_slice][60]/divisors[i], linewidth=line_width, color=colours[i], linestyle='-', label=title)
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
plt.savefig(os=path.join(args.output, f"{args.output_filename}_axial.png"))

# In[62]:

# Define consistent font size and colormap
font_size = 14
colormap = 'viridis'
coronal_slice = 55

# Set the maximum intensity for color normalization
vmax = 100#max(measured_data.max(), simind_total.max(), stir_forward_projection.max())

data_list = [
    ((simind_total).as_array(), "SIMIND total"),
    (measured_data.as_array(), "measured"),
    ((simind_true).as_array(), "SIMIND true"),
    ((stir_forward_projection).as_array(), "STIR forward"),
    ((simind_scatter).as_array(), "SIMIND scatter"),
]

# Create a figure and a GridSpec with 3 rows
fig = plt.figure(figsize=(len(data_list)*4,7*2,))
gs = GridSpec(3, len(data_list), height_ratios=[2, 0.15, 3])  # Adjusted GridSpec for clarity

# Create image subplots in the first row
ax_images = [fig.add_subplot(gs[0, i]) for i in range(len(data_list))]

for i, (data, title) in enumerate(data_list):
    im = ax_images[i].imshow(data[0, :, coronal_slice]/divisors[i], vmin=0, vmax=vmax, cmap=colormap)
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
    ax_line.plot(data[0, 60, coronal_slice]/divisors[i], linewidth=line_width, color=colours[i], linestyle='-', label=title)
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
plt.savefig(os.path.join(args.output, f"{args.output_filename}_coronal.png"))

