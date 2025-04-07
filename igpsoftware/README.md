# IGP software folder

The `igpsoftware` folder contains the Matlab code (functions) for the integrated processing.

## Contents

```
igpsoftware/
	stmmain/        > main functions
	stmutil/        > Space-Time Matrix utilities toolbox
	crsutil/        > Coordinate Reference System utilities toolbox
	mht/            > Multiple hypothesis testing toolbox (InSAR)
	proj/           > Projection toolbox
	rdnaptrans/     > RDNAP transformation toolbox (Dutch reference systems)
	tseries2/       > GNSS timeseries toolbox
	igptoolbox.cfg	
	igpimport.m
```

The function `igpimport.m` and file `igptoolbox.cfg` is used by the scripts in the `igpprojects` folder to import the above toolboxes. 
The file `igptoolbox.cfg` is used to link toolbox names in the processing scripts, which call `igpimport.m`, to a specific folder.

The toolboxes `stmmain` and `stmutil` are specific for the IGP software. The other toolboxes, `crsutil`, `tseries2`, `mht`, `proj` and `rdnaptrans`
are general purpose toolboxes that are also used by other software projects. For convenience, the code for the general purpose toolboxes are copied into
the IGP software folder, but you can have `igptoolbox.cfg` also point to other locations.

## One time action upon installation

Upon intitial installation you must modify the file `igptoolbox.cfg`. The file `igptoolbox.cfg` is used to link toolbox names in the scripts to a specific folder; 
this can be a subfolder in `igpsoftware` or any other folder. 
This facilitates easy configuration for test using modified toolboxes without additionally modifying the calling scripts.
