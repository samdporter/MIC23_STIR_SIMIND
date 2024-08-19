python STIR_demonstration.py --total_activity=258.423 --time_per_projection=43 \
--photon_multiplier=0.001 --photopeak_energy=208 --window_lower=187 --window_upper=229 \
--source_type="lu177" --collimator="G8-MEGP" --kev_per_channel=10 \
--max_energy=498 --mu_map_path="data/Lu177/registered_CTAC.hv" \
--image_path="data/Lu177/osem_reconstruction_postfilter_555.hv" \
--measured_data_path="data/Lu177/SPECTCT_NEMA_128_EM001_DS_en_1_Lu177_EM.hdr" \
--output_dir="simind_output" --input_smc_file_path="input/input.smc" \
--scoring_routine=1 --collimator_routine=1 --photon_direction=3 --crystal_thickness=7.25 \
--crystal_half_length_radius=196.8 --crystal_half_width=255.85