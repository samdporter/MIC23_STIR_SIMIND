#!/usr/bin/env python3
import argparse
import glob
import os

import matplotlib.pyplot as plt
from sirf.STIR import (
    AcquisitionData,
    SPECTUBMatrix,
    AcquisitionModelUsingMatrix,
    ImageData,
)


def get_spect_data(path):
    spect_data = {}
    spect_data["acquisition_data"] = AcquisitionData(
        os.path.join(path, "peak.hs")
    )
    spect_data["attenuation"] = ImageData(
        os.path.join(path, "umap_zoomed.hv")
    )
    # attn_arr = spect_data["attenuation"].as_array()
    # attn_arr = np.flip(attn_arr, axis=-1)
    # spect_data["attenuation"].fill(attn_arr)
    try:
        spect_data["initial_image"] = ImageData(
            os.path.join(path, "initial_image.hv")
        ).maximum(0)
    except Exception:
        spect_data["initial_image"] = ImageData(
            os.path.join(path, "template_image.hv")
        )
        spect_data["initial_image"].fill(1)

    return spect_data


def get_spect_am(spect_data, keep_all_views_in_cache=False):
    spect_am_mat = SPECTUBMatrix()
    spect_am_mat.set_attenuation_image(spect_data["attenuation"])
    spect_am_mat.set_keep_all_views_in_cache(keep_all_views_in_cache)
    spect_am_mat.set_resolution_model(0.9323, 0.03, False)
    spect_am = AcquisitionModelUsingMatrix(spect_am_mat)
    return spect_am


def main():
    parser = argparse.ArgumentParser(
        description="Compute mean scatter image from SIMIND scatter outputs."
    )
    parser.add_argument(
        '--input_dir', type=str, required=True,
        help="Directory containing scatter files"
    )
    parser.add_argument(
        '--data_dir', type=str, default=None,
        help="Directory containing data files"
    )
    parser.add_argument(
        '--scatter_pattern', type=str, default="*_sca_w1.hs",
        help="Filename pattern for scatter files (default: '*_sca_w1.hs')"
    )
    parser.add_argument(
        '--total_pattern', type=str, default="*_tot_w1.hs",
        help="Filename pattern for total files (default: '*_tot_w1.hs')"
    )
    parser.add_argument(
        '--image_pattern', type=str, default="recon_osem.hv",
        help="Filename pattern for image files (default: 'recon_osem.hv')"
    )
    parser.add_argument(
        '--output_file', type=str, required=True,
        help="Output file for mean scatter image"
    )
    parser.add_argument(
        '--delete_files', action='store_true',
        help="Delete scatter files after computing mean scatter image"
    )
    args = parser.parse_args()

    # Sum scatter files
    scatter_files = glob.glob(os.path.join(args.input_dir, args.scatter_pattern))
    if not scatter_files:
        raise ValueError(
            f"No files found in {args.input_dir} matching pattern {args.scatter_pattern}"
        )

    for i, file in enumerate(scatter_files):
        scatter = AcquisitionData(file)
        if i == 0:
            sum_scatter = scatter.get_uniform_copy(0)
        sum_scatter += scatter

    # Sum total files
    total_files = glob.glob(os.path.join(args.input_dir, args.total_pattern))
    if not total_files:
        raise ValueError(
            f"No files found in {args.input_dir} matching pattern {args.total_pattern}"
        )

    for i, file in enumerate(total_files):
        total = AcquisitionData(file)
        if i == 0:
            sum_total = total.get_uniform_copy(0)
        sum_total += total

    # Compute trues projection
    sum_trues = sum_total - sum_scatter

    spect_data = get_spect_data(args.data_dir)
    spect_am = get_spect_am(spect_data, keep_all_views_in_cache=False)
    spect_am.set_up(spect_data["acquisition_data"], spect_data["initial_image"])

    image = ImageData(os.path.join(args.input_dir, args.image_pattern))
    forward = spect_am.forward(image)

    attenuation_image = image.clone()
    forward_attenuation = spect_am.forward(attenuation_image)
    thresh = 0.01 * forward_attenuation.max()
    forward_attenuation.fill(forward_attenuation.as_array() > thresh)

    # Mask trues and forward projections
    sum_trues_masked = sum_trues.clone()
    sum_trues_masked *= forward_attenuation
    sum_trues_counts = sum_trues_masked.sum()

    forward_masked = forward.clone()
    forward_masked *= forward_attenuation
    forward_counts = forward_masked.sum()

    scatter_scaling = forward_counts / sum_trues_counts
    print(f"Scatter scaling factor: {scatter_scaling}")

    mean_scatter = sum_scatter * scatter_scaling
    mean_scatter.write(args.output_file)
    print(
        f"Mean scatter image computed from {len(scatter_files)} files and "
        f"written to {args.output_file}"
    )

    if args.delete_files:
        # Delete binary data files
        for pattern in ["*_sca_w1.a00", "*_air_w1.a00", "*_tot_w1.a00"]:
            for file in glob.glob(os.path.join(args.input_dir, pattern)):
                os.remove(file)


if __name__ == '__main__':
    main()
