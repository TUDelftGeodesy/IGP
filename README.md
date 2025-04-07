# Integrated Geodetic Processing (IGP)

**Release 1.0b1 d.d. December 2nd, 2024**

This repository contains version 1.0 of the NAM Integrated Geodetic Processing
software (IGP).

The main folder structure is

| folder       | content                                      |
| ------------ | ---------------------------------------------|
| igpsoftware  | actual Matlab software                       |
| igptemplates | processing templates (Matlab scripts)        |
| igpprojects  | processing results (empty)                   |
| igpdata      |  staging area for the input files (empty)    |

Actual processing is done within igpprojects. To start a new project in igpproject:

1. copy the raw data files to `igpdata`
2. copy one of the template folders from `igptemplates` to `igpprojects`
3. modify the copied template in `igpprojects`
4. process the project in `igpprojects`

See the IGP User Manual for more instructions.
