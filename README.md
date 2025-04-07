# Integrated Geodetic Processing (IGP)

**Release 1.0b1 d.d. December 2nd, 2024**

This repository contains version 1.0 of the NAM Integrated Geodetic Processing software (IGP) for the integrated processing of multiple geodetic datasets, acquired by levelling, GNSS and/or InSAR. 

At the core of the processing is the so-called `Space-Time Matrix (STM)` format and processing modules that operate on space-time matrices. 
All data is initially converted into the generic `STM` format. 
The `STM` datasets can be both input to further processing steps, sometimes accepting multiple inpust STM datasets, and output from processing steps. 
This allows the user to set-up a processing chain to handle a variety of geodetic data analysis and integration problems.

The processing modules are divided into

1. **Initialization**: Transforms the original datasets into the common Space-Time Matrix (STM) dataformat. 
  
2. **Selection, Decomposition and Reduction**: Reduces the GNSS CORS and InSAR data to common evaluation points and epochs. 
   
   a. Selection of evaluation points and epochs from all available input space time matrix datasets.

   > **Evaluation points** are all points present in levelling, GNSS CORS and campaign datasets, plus a densification based on the InSAR datasets.
   > **Evaluation epochs** are all epoch present present in levelling and GNSS campaing datasets, plus a densification based on the available InSAR datasets and GNSS CORS datasets.
   
   b. Decomposition and reduction.

   > The GNSS CORS and InSAR datasets are reduced to the evaluation points and epochs. For the GNSS CORS datasets the signals are also decomposed into a trend signal, harmonics, temperature and atmospheric loading effects, steps and residual signal. For the InSAR data a decomposition into shallow and deep signals is carried out.
   
3. **Integration**: Integrates multiple input STM datasets at the common evaluation points and epochs with statistical testing. The output is a single integrated space time matrix dataset.  

4. **Prediction**: Least-squares prediction to any desired location and time.

9. **Output and visualization**: The software comes with a number of utility tools to display and visualize the contents of space time matrices. These tools include a simple display of the space time datasets meta data, plotting of the space time dataset network map, time series, estimate velocities and covariance matrix, and plotting and gridding of predicted displacements.  

## Software Installation

Clone this repository at a folder of your choice. This results in the following four directories:

| folder       | content                                      |
| ------------ | ---------------------------------------------|
| igpsoftware  | actual Matlab software                       |
| igptemplates | processing templates (Matlab scripts)        |
| igpprojects  | processing results (empty)                   |
| igpdata      | staging area for the input files (empty)     |

To complete the installation, 

a. Modify the file `igpinit.m` in the `igpprojects` folder. The Matlab path in the file `igpinit.m`  must be set to the main `igpsoftware` folder, and the structure variable `globalAttributesDefaults` must be defined with default values for the `globalAttributes`.

b. Modify the file `igptoolbox.cfg` in the `igpsoftware` folder. The file `igptoolbox.cfg` is used to link toolbox names in the scripts to a specific folder; this can be a subfolder in `igpsoftware` or any other folder. This facilitates easy configuration for test using modified toolboxes without additionally modifying the calling scripts.

In principle, the IGP software can be run on any platform. However, especially the initialization and reduction step for InSAR data require a considerable amount of RAM. Therefore, use of a high-performance system may be needed (for these steps). Furthermore, the use of relatively new in-build Matlab functions and toolboxes is avoided as much as possible, so the software should also run with older versions of Matlab

 
## Starting a new project

Actual processing is done within `igpprojects`. To start a new project in `igpprojects`:

1. copy the raw data files to `igpdata` and copy one of the template folders from `igptemplates` to `igpprojects`
2. modify the copied template in `igpprojects`
3. process the project in `igpprojects`

The difference between `igptemplates` and `igpprojects` is that `igptemplates` contains initial example Matlab scripts,  whereas the subfolders in `igpprojects` are working directories which contain in addition to the (modified) Matlab scripts input, intermediate, output and log files. To start a new project in `igpprojects`, first copy a template from `igptemplates`, modify the scripts and then run.

More instructions can be found in the [IGP User Manual](IGP_user_manual_v1-0.pdf).

