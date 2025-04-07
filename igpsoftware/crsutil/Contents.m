% CRSUTIL Coordinate and Time Reference System Toolbox.
% Version 1.3 (9 October 2017).
%
% Coordinate transformations (ECEF reference frame)
%
%   xyz2plh     - Cartesian Coordinates to Ellipsoidal coordinates
%   plh2xyz     - Ellipsoidal coordinates to Cartesian Coordinates
%   inqell      - Semi-major axis, flattening and GM for various ellipsoids
%   plh2str     - Phi,lambda,h to a string with dd mm ss.ssss notation
%   str2plh     - Reads Phi, lambda and h from a string
%   xyz2neu     - North, East, Up (dN, dE, dU) to Cartesian delta's (dX, dY, dZ)
%   neu2xyz     - Cartesian delta's (dX, dY, dZ) to North, East, Up (dN, dE, dU)
%   plh2neu     - Ellipsoidal (Lat,Lon,Hgt) to North,East,Up (dN, dE, dU)
%   xyz2zas     - Cartesian coordinates to Zenith angle, azimuth and distance
%   zas2xyz     - Zenith angle, azimuth and distance to cartesian coordinates
%
% Time conversions (with Matlab date numbers)
%
%   date2mjd    - Matlab date to Modified Julian Date
%   mjd2date    - Modified Julian Date to Matlab date number
%   date2gps    - Matlab date to GPS second of week and GPS week
%   gps2date    - GPS second of week and GPS week to Matlab date number
%   date2dyear  - Matlab date to decimal year
%   dyear2date  - decimal year to Matlab date number
%
% Time conversions (other/legacy functions)
%
%   ymd2mjd     - [year month day] to Modified_Julian_Date
%   mjd2ymd     - Modified_Julian_Date to [year month day]
%   gps2mjd     - [GPSweek GPSsecond_in_week] to Modified_Julian_Date
%   mjd2gps     - Modified_Julian_Date to [GPSweek GPSsecond_in_week]
%   ymd2dyear   - [year month day hour min sec] to decimal year
%   dyear2ymd   - decimal year to [year month day hour min sec]
%   hms2sec     - [hour  minute second] to second_in_day
%   sec2hms     - second_in_day to [hour minute second]
%   hms2str     - [hour  minute second] to string with 'hh:mm:ss.sss'
%   str2hms     - string with 'hh:mm:ss.sss' to [hour  minute second]
%   gps2str     - GPSsecond_in_week to string with 'dd-mon-yy hh:mm:ss.sss'
%   str2gps     - string with 'dd-mmm-yy hh:mm:ss.sss' to GPSsecond_in_week
%   sec2str     - Second_in_day to string with 'hh:mm:ss.sss'
%   str2sec     - string with 'hh:mm:ss.sss' to second_in_day
%   ymd2str     - [year month day ] to string with 'dd-mmm-yy'
%   str2ymd     - string with 'dd-mmm-yy' or 'dd-mm-yy' to [year month day]
%   gpsdate     - prints the GPS week number, second in the week and day of year
%
% UT1 to GMST, and ECI/ECEF, conversions
%
%   ut2gmst    - Compute Greenwich Mean Siderial Time from UT1
%   ecef2eci   - Convert position and velocity from ECEF to ECI reference frame
%   eci2ecef   - Convert position and velocity from ECI to ECEF reference frame
%
% Keplerian elements
%
%   vec2orb     - Convert inertial state vector into Keplerian elements
%   orb2vec     - Convert Keplerian elements into inertial state vector
%   kepler      - Compute mean anomaly from eccentric anomaly (Kepler's equation)
%   keplernu    - Compute mean anomaly from true anomaly (Kepler's equation)
%   keplerm     - Compute eccentric/true from mean anomaly solving Kepler's eqn
%
% Miscellaneous
%
%   prtcrd      - Print table with coordinates and optional co-variances
%   covreformat - Reformat a co-variance matrix
%
% (c) Hans van der Marel, Delft University of Technology, 1995-2017.