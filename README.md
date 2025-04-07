# Integrated Geodetic Processing (IGP)

**Release 1.0b1 d.d. December 2nd, 2024**

This repository contains version 1.0 of the NAM Integrated Geodetic Processing software (IGP). 

The IGP software enables an integrated processing of various geodetic datasets, acquired by levelling, GNSS and/or InSAR. 

At the core of the processing is the so-called `Space-Time Matrix (STM)` format. All data is initially converted into the generic `STM` format. The `STM` datasets can be both input to further processing steps, sometimes accepting multiple inpust STM datasets, and output from processing steps. This allows the user to set-up a processing chain to handle a variety of geodetic data analysis and integration problems.

The processing modules can divided into

1. **Initialization**: Transforms the original datasets into the common Space-Time Matrix (STM) dataformat. 
  
2. **Reduction**: Reduces the GNSS CORS and InSAR data to common evaluation points and epochs. This consists of two steps
   
   a. Selection of evaluation points and epochs from all available input space time matrix datasets.

      The output is a *project* space time matrix containing the selected points and epochs. **Evaluation points** are all points present in levelling, GNSS CORS and campaign datasets, plus a densification based on the InSAR datasets. **Evaluation epochs** are all epoch present present in levelling and GNSS campaing datasets, plus a densification based on the available InSAR datasets and GNSS CORS datasets.
   
   b. Decomposition and reduction.

   In this step the GNSS CORS and InSAR datasets are reduced to the evaluation points and epochs. For the GNSS CORS datasets the signals are first decomposed into a trend signal, harmonics, temperature and atmospheric loading effects, steps and residual signal. The reduced dataset contains typically only the trend and residual signal. For the InSAR data a decomposition into shallow and deep signals is carried out.
   
4. **Integration**: Integrates multiple input STM datasets at the common evaluation points and epochs with statistical testing. The output is a single integrated space time matrix dataset.  The Integration step can be performed on a certain Region and/or Period of Interest.

5. **Prediction**:

The software comes with a number of utility tools to display and visualize the contents of space time matrices. These tools include a simple display of the space time datasets meta data, plotting of the space time dataset network map, time series, estimate velocities and covariance matrix, and plotting and gridding of predicted displacements.  

## Software Installation

At the location of choice, clone this repository. This results in the following four directories:

| folder       | content                                      |
| ------------ | ---------------------------------------------|
| igpsoftware  | actual Matlab software                       |
| igptemplates | processing templates (Matlab scripts)        |
| igpprojects  | processing results (empty)                   |
| igpdata      | staging area for the input files (empty)    |

To complete the installation, you must modify the file `igpinit.m` in the `igpprojects` directory.
In the file `igpinit.m` the Matlab path must be set to the main `igpsoftware` location, and the structure variable `globalAttributesDefaults` must be defined with default values for the `globalAttributes`.

In principle, the IGP software can be run on any platform. However, especially the initialization and reduction step for InSAR data require a considerable amount of RAM. Therefore, use of a high-performance system may be needed (for these steps). Furthermore, the use of relatively new in-build Matlab functions and toolboxes is avoided as much as possible, so the software should also run with older versions of Matlab

 
## Processing

Actual processing is done within igpprojects. To start a new project in igpproject:

1. copy the raw data files to `igpdata`
2. copy one of the template folders from `igptemplates` to `igpprojects`
3. modify the copied template in `igpprojects`
4. process the project in `igpprojects`

See the IGP User Manual for more instructions.

