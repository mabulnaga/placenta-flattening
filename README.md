# placenta-flattening
A MATLAB algorithm for volumetric mesh parameterization. Developed for mapping a placenta segmentation derived from an MRI image to a flattened template for visualization. The code can work on NIFTI images or MATLAB matrices containing imaging information.

### Requirements
- MATLAB
- GPU with CUDA compute capability
- MATLAB Packages:
    - [iso2mesh](http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download): mesh generating toolbox
    - [geometry processing toolbox](https://github.com/alecjacobson/gptoolbox)
    - [medical imaging MATLAB toolbox](https://github.com/adalca/matlib) (only needed if working with NIFTI files)
    
Add the MATLAB packages to the working path.

### Usage
main(grayImage, segImage): input a grayscale MRI image and the corresponding binary segmentation, where voxels labeled '1' correspond to the placenta. The inputs grayImage, segImage can either be full path locations of NIFTI image files, or image matrices. The grayImage input can be a 3D MRI volume, or a 4D series of MRI volumes. The script outputs the flattened meshes and images containing the mapped intensities.

### Development
Please contact Mazdak Abulnaga, abulnaga@mit.edu.

### Citing and Paper

Please consider citing our [paper](https://arxiv.org/pdf/1903.05044.pdf)
```
@inproceedings{abulnaga2019placenta,
title={Placental Flattening via Volumetric Parameterization},
author={Abulnaga, S. Mazdak and Abaci Turk, Esra and Bessmeltsev, Mikhail and Grant, P. Ellen and Solomon, Justin and Golland, Polina},
booktitle={Medical Image Computing and Computer Assisted Intervention -- MICCAI 2019},
year={2019},
pages={39--47},
}
```
