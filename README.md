# placenta-flattening
A MATLAB algorithm for volumetric mesh parameterization. Developed for mapping a placenta segmentation derived from an MRI image to a flattened template for visualization. The code can work on NIFTI images or MATLAB matrices containing imaging information.

![alt text](https://github.com/mabulnaga/placenta-flattening/blob/master/flattening_flowchart.png)

### Requirements
- MATLAB [![View placenta-flattening on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/114395-placenta-flattening)
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
If you use this method or some parts of the code, please consider citing one of our papers. 

Our journal paper develops additional template models and provides extensions to improve robustness, an expanded evaluation on a significantly larger dataset, and experiments demonstrating utility for clinical research. [eprint arXiV:2111.07900](https://arxiv.org/abs/2111.07900)
```
@ARTICLE{abulnaga2021placenta,
  author={Abulnaga, S. Mazdak and Abaci Turk, Esra and Bessmeltsev, Mikhail and Grant, P. Ellen and Solomon, Justin and Golland, Polina},
  journal={IEEE Transactions on Medical Imaging}, 
  title={Volumetric Parameterization of the Placenta to a Flattened Template}, 
  year={2021},
  volume={},
  number={},
  pages={1-1},
  doi={10.1109/TMI.2021.3128743}}
```

The MICCAI conference paper develops the parallel planes template and validates on a smaller dataset. [eprint arXiV:1903.05044](https://arxiv.org/abs/1903.05044)
```
@inproceedings{abulnaga2019placenta,
title={Placental Flattening via Volumetric Parameterization},
author={Abulnaga, S. Mazdak and Abaci Turk, Esra and Bessmeltsev, Mikhail and Grant, P. Ellen and Solomon, Justin and Golland, Polina},
booktitle={Medical Image Computing and Computer Assisted Intervention -- MICCAI 2019},
year={2019},
pages={39--47},
}
```
