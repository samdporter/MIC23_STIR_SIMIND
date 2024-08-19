python STIR_demonstration.py --total_activity=187 --time_per_projection=20 \
--photon_multiplier=100 --photopeak_energy=150 --window_lower=75 --window_upper=225 \
--source_type="y90" --collimator="ma-megp" --kev_per_channel=25 \
--max_energy=2180 --mu_map_path="/home/sam/data/phantom_data/SPECT/umap_zoomed.hv" \
--image_path="/home/sam/data/phantom_data/SPECT/ellipsoid_image_s.hv" \
--measured_data_path="/home/sam/data/phantom_data/SPECT/peak_1_projdata__f1g1d0b0.hs" \
--output_dir="../simind_output" --input_smc_file_path="../input/input.smc" \
--scoring_routine=1 --collimator_routine=1 --photon_direction=-90 --crystal_thickness=15.9 \
--crystal_half_length_radius=235 --crystal_half_width=285
