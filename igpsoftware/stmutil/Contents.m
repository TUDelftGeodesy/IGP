% STM  Matlab Toolbox STMUTIL
%
% The STMUTIL toolbox contains high-level functions for reading, writing, 
% query, process, print and plot Space Time Matrix datasets.
%
% Top-level (command line) printing and plotting of space time matrices:
% 
%    stmdisp           - Display basic information on space time matrix datasets 
%    stmdiff           - Diff (compare) two space time matrix datasets
%
%    stmplotmap        - Plot map with the points in the stm dataset 
%    stmplotmapbyepoch - Plot map of the space time matrix points by epoch.
%    stmplotprojectmap - Plot map of the evaluation points 
%    stmplotseries     - Plot timeseries from space time matrix datasets 
%    stmplotvel        - Plot velocities estimated by stmvelocity 
%    stmplotcov        - Plot space time matrix covariance matrix
%    stmresplot        - Plot residuals and test statistics from stmintegrate
%
% The top-level functions operate on space time matrices directly and are typically
% used to print and plot information from space time matrix datasets from the
% command line or scripts. 
%
% High-level creation, reading and writing of space time matrices:
%
%    stm               - Create a Space Time Matrix dataset structure or object
%    stmread           - Read Space Time Matrix dataset from file
%    stmwrite          - Write Space Time Matrix structure to file
%
% High-level Space time matrix operations:
%
%    stminterpolate    - Interpolate space time matrix to new set of points and epochs  
%    stmresiduals      - Residual Testing for Integrated processing
%    stmvelocity       - Velocity estimation for space time matrix
%    stmprediction     - Spatio-temporal prediction based on space time matrix
%    stmrestore        - Restore of estimated trend (velocity) and predicted signal
%    stmclustering     - Clustering of points in space time matrix
%
% High-level Stochastic model definition, interogation and handling:
%
%    stmstochmodel     - Create covariance matrix from STM model specification  
%    printStochModel   - Print stochastic model for space time datasets
%    updStochModel     - Update stochastic model for space time datasets
%
% Low-level functions:
%
%    stmcheckarguments - Check input arguments for the main modules
%    stmprintopt       - Pretty print STM options
%    getinputfilenames - Read input files names from file.
%    
%    roi2poly          - Convert Region Of Interest (ROI) file/string/numerical into a polygon.
%    getepochmask      - Compute epoch mask for Period Of Interest (POI).
%    getpntmask        - Compute point mask from Region Of Interest (ROI) file/string/numerical.
%
%    olcencode         - Encode Open Location Code (plus code)
%    plh2neusp         - Ellipsoidal (Lat,Lon,Hgt) to local coordinates (North,East,Up)
%    neu2plhsp         - Local coordinates (North,East,Up) to Ellipsoidal (Lat,Lon,Hgt)
%
%    chol2stmcov       - Compute covariance matrix from Cholesky factor
%    symdiags          - Store covariance matrix in diagonals format
%
%    stmtrans          - Compute Space Time Matrix S-Transformation basis
%    stmtransapply     - Apply Space Time Matrix S-Transformation
%
%    crstrans          - Coordinate reference system transformations.
%
% The STM dataset is stored as a Matlab V7.3 mat file. For details on the dataset 
% stucture see the help of the STM function.
%
% (c) Hans van der Marel, Freek van Leijen, Delft University of Technology, 2020-2024.
