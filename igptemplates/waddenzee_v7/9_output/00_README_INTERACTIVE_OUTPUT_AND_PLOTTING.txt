Interactive output and plotting
-------------------------------

Stmutil contains several top-level functions for printing and plotting of space 
time matrices:
 
    stmdisp           - Display basic information on space time matrix datasets 
    stmdiff           - Diff (compare) two space time matrix datasets

    stmplotmap        - Plot map with the points in the stm dataset 
    stmplotprojectmap - Plot map of the evaluation points 
    stmplotseries     - Plot timeseries from space time matrix datasets 
    stmplotcov        - Plot space time matrix covariance matrix
    stmplotvel        - Plot velocities estimated by stmvelocity 
    stmresplot        - Plot residuals and test statistics from stmintegrate

The top-level functions operate on space time matrices directly and are typically
used to print and plot information from space time matrix datasets from the
command line or scripts. 

To prepare plot backgrounds run once the script (in this folder):
 
    igp_background    - Creates igp-background.mat file from shapefiles in stmutil/basemaps
    igp_unstablearea  - Creates igp-unstablearea.mat file from unstablearea.coo

For examples of the top-level plotting functions see the script

    stmplot_examples  - Examples for plotting of stm files 
  
  