### These can now be found in SIRF utils


import os
import sirf.STIR as spect
import numpy as np

def create_stir_image(matrix_dim:list, voxel_size:list):
    '''
    Creates a uniform (zeros) STIR ImageData object given specified parameters.

    Parameters:
    matrix_dim (list of int): A three element list containing the matrix size for each dimension of the image.
    voxel_size (list of float): A three element list describing the voxel size in the image (mm).

    Returns [ImageData]: The ImageData object.

    '''
    img = np.zeros(matrix_dim, dtype=np.float32)

    header = {}
    header['!INTERFILE'] = ''
    header['!imaging modality'] = 'nucmed'
    header['!version of keys'] = 'STIR4.0'

    header['!GENERAL DATA'] = ''
    header['!name of data file'] = 'temp.v'

    header['!GENERAL IMAGE DATA'] = ''
    header['!type of data'] = 'Tomographic'
    header['imagedata byte order'] = 'LITTLEENDIAN'

    header['!SPECT STUDY (general)'] = ''
    header['!process status'] = 'reconstructed'
    header['!number format'] = 'float'
    header['!number of bytes per pixel'] = '4'
    header['number of dimensions'] = str(np.size(img.shape))
    header['matrix axis label [1]'] = 'x'
    header['matrix axis label [2]'] = 'y'
    header['matrix axis label [3]'] = 'z'
    header['!matrix size [1]'] = str(matrix_dim[2])
    header['!matrix size [2]'] = str(matrix_dim[1])
    header['!matrix size [3]'] = str(matrix_dim[0])
    header['scaling factor (mm/pixel) [1]'] = str(voxel_size[2])
    header['scaling factor (mm/pixel) [2]'] = str(voxel_size[1])
    header['scaling factor (mm/pixel) [3]'] = str(voxel_size[0])
    header['number of time frames'] = '1'

    header['!END OF INTERFILE'] = ''

    line = 0
    header_path = os.path.join('temp.hv')
    with open(header_path, 'w') as f:
        for k in header.keys():
            if k.islower() or line == 0:
                tempStr = str(k)+' := '+str(header[str(k)])+'\n'
                line +=1
            else:
                tempStr = '\n'+str(k)+' := '+str(header[str(k)])+'\n'
            f.write(tempStr)
            #print(k, ":=", header[str(k)])

    f.close()

    raw_file_path = os.path.join('temp.v')
    img.tofile(raw_file_path)

    print('Image written to: ' + header_path)

    template_image = spect.ImageData(header_path)
    os.remove(header_path)
    os.remove(raw_file_path)

    return template_image


def create_stir_acqdata(proj_matrix:list, num_projections:int, pixel_size:list):
    '''
    Creates a uniform (zeros) STIR AcquisitionData object given specified parameters.

    Parameters:
    proj_matrix (list of int): A two element list containing the matrix size for each dimension of the projections.
    num_projections (int): The number of projections in the acquisition data file.
    pixel_size (list of float): A two element list describing the pixel size in the projections (mm).

    Returns [AcquisiitonData]: The AcquisitionData object.

    '''
    acq = np.zeros((1, proj_matrix[0], num_projections, proj_matrix[1]), dtype=np.float32)

    header = {}
    header['!INTERFILE'] = ''
    header['!imaging modality'] = 'NM'
    header['name of data file'] = 'temp.s'
    header['!version of keys'] = '3.3'

    header['!GENERAL DATA'] = ''

    header['!GENERAL IMAGE DATA'] = ''
    header['!type of data'] = 'Tomographic'
    header['imagedata byte order'] = 'LITTLEENDIAN'

    header['!SPECT STUDY (General)'] = ''
    header['!number format'] = 'float'
    header['!number of bytes per pixel'] = '4'
    header['!number of projections'] = str(num_projections)
    header['!extent of rotation'] = '360'
    header['process status'] = 'acquired'

    header['!SPECT STUDY (acquired data)'] = ''
    header['!direction of rotation'] = 'CW'
    header['start angle'] = '180'
    header['orbit'] = 'Circular'
    header['Radius'] = '200'

    header['!matrix size [1]'] = str(proj_matrix[0])
    header['scaling factor (mm/pixel) [1]'] = str(pixel_size[0])
    header['!matrix size [2]'] = str(proj_matrix[1])
    header['scaling factor (mm/pixel) [2]'] = str(pixel_size[1])

    header['!END OF INTERFILE'] = ''

    header_path = os.path.join('temp.hs')
    with open(header_path, 'w') as f:
        for k in header.keys():
            tempStr = str(k)+' := '+str(header[str(k)])+'\n'
            f.write(tempStr)

    f.close()

    raw_file_path = os.path.join('temp.s')
    acq.tofile(raw_file_path)

    print('Acquisition Data written to: ' + header_path)

    template_acqdata = spect.AcquisitionData(header_path)
    os.remove(header_path)
    os.remove(raw_file_path)

    return template_acqdata