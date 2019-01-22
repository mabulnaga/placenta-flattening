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
main(grayImage, segImage): input a grayscale MRI image and the corresponding binary segmentation, where voxels labeled '1' correspond to the placenta. The inputs grayImage, segImage can either be full path locations of NIFTI image files, or image matrices. The script outputs the flattened meshes and images containing the mapped intensities.

### Development
Please contact Mazdak Abulnaga, abulnaga@mit.edu.

### Paper
In progress.
