python STIR_demonstration.py --initial_activity=177 --time_per_projection=32 \
--photon_multiplier=1 --photopeak_energy=140.5 --window_lower=126.459 --window_upper=154.561 \
--source_type="tc99" --collimator="ma-lehr" --kev_per_channel=1\
--max_energy=160 --mu_map_path="../data/Tc99m/umap_resampled.hv" \
--image_path="../data/Tc99m/osem_reconstruction_postfilter_555.hv" \
--measured_data_path="../data/Tc99m/peak_stir_en_1_Primary.hdr" \
--output_dir="../simind_output" --input_smc_file_path="../input/input.smc" \
--scoring_routine=1 --collimator_routine=0 --photon_direction=2 --crystal_thickness=15.9 \
--crystal_half_length_radius=235 --crystal_half_width=285 \
--half_life=6