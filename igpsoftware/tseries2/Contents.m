% TSERIES GNSS Time Series Analysis Toolbox.
% Version 1.2b (4 March 2019). Updated 19 September 2019, 10 March 2021.
%
% Time series fitting and analysis
%
%   tseriesfit            - Fit time series components.
%   tseriesanalysis       - Analysis of timeseries of GPS station coordinates.
%   tseriesperiodogram    - Plot timeseries periodogram for one or more stations.
%   tseriesresidualstack  - Compute and plot residual stack from fitted timeseries.
%   tseriescmfit          - Compute the common mode fit for GPS timeseries.
%   tseriescmeval         - Evaluate common mode fit.
%   tseriescmstack        - Converts residual stack and common mode into a timeseries.
%
% Plotting and printing
%
%   tseriesplot           - Plot timeseries components.
%   tseriesplotcomp       - Plot timeseries components for multiple stations, new version, 
%                           not backwards compatible with tseriesplotcomponent
%   tseriesplotcomponent  - Plot timeseries components for multiple stations.
%   tseriesplotcomponent2 - Plot timeseries components for two timeseries and for multiple stations.
%   tseriesplotmap        - Plot map with representation of the estimated parameters.
%   tseriessummary        - Print summary of estimated parameters and standard deviations.
%
% Time series output
%
%   tseriestxtwrite       - Write timeseries to an ascii text file.
%   tseriesxlswrite       - Write timeseries to an excel file.
%
% Meteo data support
%
%   getmeteo              - Read temperature and pressure data from KNMI hourly meteo files.
%
% Event management
%
%   getEvents             - Get events from csv text file.
%   prtEvents             - Print the list of events for an event or time series structure.
%   prtSteps              - Print list of steps with estimated step size.
%   trimEvents            - Trim the list of events for a given data period.
%   addEvent              - Add event to list of events from a time series structure.
%
% Filtering, interpolation and outlier/step detection
%
%   tseriesoutlier
%   rmafilt               - Robust moving average filter with outlier and step detection
%   mafilt1               - Moving Average interpolation
%   mmfilt1               - Moving Median interpolation
%
% Miscellaneous
%
%   fig2subplot           - Copies figures to subplots in new figure. 
%   savepng               - Save figure(s) as png file.
%
%   pltarrow              - plot an arrow
%   pltellips             - plot an error ellipse
%   pltseries1            - Elementary timeseries plot with annotations
%   pltseries1b           - Elementary timeseries plot with annotations
%   pltspectrum           - Plot GPS Lomb-Scargle Periodogram
%
%   plomb                 - Lomb-Scargle periodogram
%
% (c) Hans van der Marel, Delft University of Technology, 2004-2021.