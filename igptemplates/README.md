# IGP templates

Within the `igpprojects` and `igptemplates` folder different **projects** can be setup. 


## Contents

A typical layout for the `igpprojects` and `igptemplates` directory is

```
igpprojects/
    groningen/
    simtests/
    waddenzee/
        0_import_gnss/
        0_import_insar/
        0_import_levelling/
        1_decompose_gnss/
        2_reduce/
        3_integrate/
        4_predict/
        9_output
    igpinit.m
```

The difference between `igptemplates` and `igpprojects` is that `igptemplates` contains initial example Matlab scripts,  whereas the subfolders in 
`igpprojects` are working directories which contain in addition to the (modified) Matlab scripts input, intermediate, output and log files. 

## Projects

In each project the multiple processing steps can be run. However, it is not necessary to run every step in single project; 
a project may also make use of intermediate outputs of other projects. For example, a project may use data imported in another project. 

## Starting a new project

To start a new project in `igpprojects`, first copy a template from `igptemplates`, modify the scripts and then run.

The folder structure in `igpprojects` is quite flexible and allows for various choices. 
For instance, if you decide to investigate an alternative processing, you can create a new project of you can create a new subfolder within a project 
(e.g. 3_integrate_alternative), or just create a new output file within an existing subfolder. 

The example layout contains a folder `9_output`. This folder contains initially a few mat files with the unstable area and coastlines for plotting. 
This folder is also intended for interactive plotting.

## One time action on installation

The `igpprojects` directory contains a single initialization file, `igpinit.m`,  with the path to the software directory (`igpsoftware`) and initial values 
for the `global attributes`. This is the file you have to modify at install time.
