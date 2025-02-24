from sirf.STIR import *
import numpy as np
import os
import matplotlib.pyplot as pl

import argparse

parser = argparse.ArgumentParser(description='Reconstruct with OSEM')

parser.add_argument('--data_path', type=str, default="/home/storage/copied_data/data/phantom_data/for_cluster/SPECT", help='data path')
parser.add_argument('--num_subsets', type=int, default=12, help='number of subsets')
parser.add_argument('--num_epochs', type=int, default=10, help='number of epochs')
# default additive path to None but expect string
parser.add_argument('--additive_path', type=str, default=None, help='additive path')
parser.add_argument('--smoothing', type=bool, default=True, help='smoothing')
parser.add_argument('--index', type=int, default=0, help='index')

def get_spect_data(path):

    spect_data = {}
    spect_data["acquisition_data"] = AcquisitionData(os.path.join(path,  "peak.hs"))
    try:
        spect_data["attenuation"] = ImageData(os.path.join(path,  "umap.hv"))
    except:
        spect_data["attenuation"] = ImageData(os.path.join(path,  "umap_zoomed.hv"))
    #attn_arr = spect_data["attenuation"].as_array()
    #attn_arr = np.flip(attn_arr, axis=-1)
    #spect_data["attenuation"].fill(attn_arr)
    spect_data["initial_image"] = ImageData(os.path.join(path,  "initial_image.hv")).maximum(0)

    return spect_data

def get_spect_am(spect_data, keep_all_views_in_cache=True):
    spect_am_mat = SPECTUBMatrix()
    spect_am_mat.set_attenuation_image(spect_data["attenuation"])
    spect_am_mat.set_keep_all_views_in_cache(keep_all_views_in_cache)
    spect_am_mat.set_resolution_model(0.9323, 0.03, False) 
    spect_am = AcquisitionModelUsingMatrix(spect_am_mat)
    if spect_data["additive"] is not None:
        spect_am.set_additive_term(spect_data["additive"])
    return spect_am

def get_reconstructor(data, acq_model, initial_image, num_subsets, num_epochs):
    recon = OSMAPOSLReconstructor()
    recon.set_objective_function(make_Poisson_loglikelihood(acq_data = data, acq_model = acq_model))
    recon.set_num_subsets(num_subsets)
    recon.set_num_subiterations(num_subsets * num_epochs)
    recon.set_up(initial_image)
    return recon

def main(data_path):

    spect_data = get_spect_data(data_path)
    if args.additive_path is not None:
        spect_data['additive'] = AcquisitionData(args.additive_path)
    else:
        spect_data['additive'] = None
    spect_am = get_spect_am(spect_data, True)
    spect_init = spect_data["initial_image"]
    spect_recon = get_reconstructor(spect_data["acquisition_data"], spect_am, spect_init, args.num_subsets, args.num_epochs)
    spect_recon.reconstruct(spect_init)

    recon_image = spect_recon.get_current_estimate()

    if args.smoothing:
        gauss = SeparableGaussianImageFilter()
        gauss.set_fwhms((5, 5, 5))
        gauss.apply(recon_image)

    return recon_image

if __name__ == "__main__":
    
    msg = MessageRedirector()
    
    args = parser.parse_args()
    spect = main(args.data_path)
    suffix = f"osem_i{args.num_epochs}_s{args.num_subsets}"
    if args.smoothing:
        suffix += "_smoothed"
    spect.write(os.path.join(args.data_path, f"recon_{suffix}_{args.index}.hv"))

