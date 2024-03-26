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
```[startVolume, mappedVolume, mappedImage] = main(grayImage, segImage)```

- grayImage: grayscale MRI image. Input can be a 3D MRI volume, or a 4D series of MRI volumes. 

- segImage: 3D binary segmentation image, where voxels labeled '1' correspond to the placenta.

Either input can be a full path location pointing to the NIFTI image files, or matrices. 

The script outputs the source mesh (startVolume), the flattened mesh (mappedVolume), and an MRI image containing the mapped intensities (mappedImage).

If you have multiple sources of MRI data corresponding to the same placenta segmentation, you can map each of these individually to the flattened space by running the command:

``` [mappedImage] = map_MRI_intensity(startVolume, mappedVolume, mriImage) ```

mriImage is a 3D or 4D matrix containing MRI signals to be mapped to the flattened space.

### Development
Please contact Mazdak Abulnaga, abulnaga@mit.edu.

### Citing and Paper
If you use this method or some parts of the code, please consider citing one of our papers. 

Our journal paper develops additional template models and provides extensions to improve robustness, an expanded evaluation on a significantly larger dataset, and experiments demonstrating utility for clinical research. [eprint arXiV:2111.07900](https://arxiv.org/abs/2111.07900)
```
@ARTICLE{abulnaga2022placenta,
  author={Abulnaga, S. Mazdak and Abaci Turk, Esra and Bessmeltsev, Mikhail and Grant, P. Ellen and Solomon, Justin and Golland, Polina},
  journal={IEEE Transactions on Medical Imaging}, 
  title={Volumetric Parameterization of the Placenta to a Flattened Template}, 
  year={2022},
  volume={41},
  number={4},
  pages={925-936},
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
