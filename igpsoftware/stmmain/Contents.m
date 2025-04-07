% STM  Matlab Toolbox STMMAIN
%
% STMMAIN contains the main building blocks for Integrated Geodetic Geodetic
% Processing (IGP) using the Space Time Matrix (STM) concept.
%
% Space Time Matrix creation and conversion:
%
%    cupido2stm        - CUPiDO v1 to Space Time Matrix format conversion.
%    gnss2stm          - Import gnss data files into space time matrix dataset
%    insar2stm         - Import InSAR data files into space time matrix dataset
%    simulate2stm      - Simulate 2-D deformation field with geodetic observations.
%
% Space Time Matrix reduction and decomposition:
%
%    stmdecomposegnss  - Decompose GNSS time series.
%
%    stmselect         - Select evaluation epochs and points for integrated processing.
%
%    stmreducecampaign - Reduce campaign datasets.
%    stmreducegnss     - Reduce GNSS time series.
%    stmreduceinsar    - Data reduction InSAR (including decomposition).
%
% Integrated processing and prediction:
%
%    stmintegrate      - Integrated processing of geodetic displacement observations.
%    stmpredict        - Prediction based on Space Time Matrix.
%
% Export of STM to alternative formats:
%
%    stmexport         - Export of STM to .csv, .png, .shp, .tiff
%
% The functions in STMMAIN have all as input and/or output one, two or 
% more space time matrices. They are typically called from scripts in
% the IGPPROJECTS subdirectories, using the syntax
%
%    FUNCTIONNAME(INPUTFILENAME(S),OUTPUTFILENAME,OPTIONS)
%
% See also STMUTIL.
%
% (c) Hans van der Marel, Freek van Leijen, Delft University of Technology, 2020-2023.
