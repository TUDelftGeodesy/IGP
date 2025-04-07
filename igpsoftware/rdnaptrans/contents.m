% RDNAPTRANS (TM) 2008 (2004) Matlab Toolbox
% Version 1.2 (12 June 2013, updated 30 April 2019)
%
% RDNAPTRANS is the official transformation between the European Terrestial
% Reference System (ETRS89), and the Dutch Rijksdriehoeks (RD) and Normal 
% Amsterdams Peil (NAP) systems.
%
% The transformation includes a 3-D similarity transformation from the ETRS89 
% datum to the Dutch datum, Stereographic Double projection onto the Bessel 
% ellipsoid in the Dutch datum, a correction grid for the RD coordinates, and 
% the NLGEOID2004 geoid. The differences do not exceed 1 cm in horizontal 
% position and 2 cm in the vertical. Reference: "A. de Bruijne, J. van Buren, 
% A. Kösters and H. van der Marel, De geodetische referentiestelsels van 
% Nederland; Definitie en vastlegging van ETRS89, RD en NAP en hun onderlinge 
% relatie, NGC, Delft, 2005". 
%
% These function provide both the current RDNAPTRANS2008 as well as the legacy 
% RDNAPTRANS2004 procedure. 
%
% High-level functions:
%
% etrs2rdnap  - Convert cartesian ETRS89 coordinates into RDNAP.
% rdnap2etrs  - Convert RDNAP coordinates into cartesian ETRS89 coordinates
% nap2etrs    - Convert NAP height into ellipsoidal ETRS89 height
% etrs2nap    - Convert ellipsoidal ETRS89 height into NAP height
% rd2etrs     - Convert RD grid coordinates into ETRS89 latitude and longitude
% etrs2rd     - Convert ETRS89 latitude and longitude into RD grid coordinates%
%
%
% Low-level functions (should not be called by the user directly):
%
% 7 Parameter similarity transformation:
%
% simtrans    - 3D-simularity transformation for cartesian coordinates
% etrs2bessel - Transformation parameters for ETRS89 to NL-Bessel datum  
% bessel2etrs - Transformation parameters for NL-Bessel to ETRS89 datum  
%
% Conversion between cartesian and geographic coordinates
%
% plh2xyz     - Ellipsoidal coordinates to Cartesian Coordinates 
% xyz2plh     - Cartesian Coordinates to Ellipsoidal coordinates
% inqell      - Semi-major axis, flattening and GM for various ellipsoids
%
% Stereographic double projection and parameters for the Netherlands
%
% stereodp    - Stereographic double projection 
% nlbessel    - Parameters for the Dutch Stereographic double projection 
%
% Grid corrections and geoid
%
% readgrd     - Read binary grid files in the format of the program Surfer(R) 
% grdint      - Grid interpolation using Catmull-Rom-Overhauser splines
%
% nlgeo04.grd - Data file with NLGEOID2004
% x2c.grd     - Data file with corrections for the x-RD coordinates
% y2c.grd     - Data file with corrections for the x-RD coordinates
%
% Test functions
%
% testrdnap   - Matlab function to test the implementation
% testgrdint  - Matlab function to test the grid-interpolation
%
%
% The Matlab functions in this package are free software: you can redistribute 
% it and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation <http://www.gnu.org/licenses/>.
%
% The data files nlgeo04.grd, x2c.grd and y2c.grd have been compiled by
% Kadaster (Apeldoorn) and Rijkswaterstaat (Delft), and are provided as is
% without modification. These files are property of Kadaster/RWS and may not
% be changed.
%
% This software is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
% (c) Hans van der Marel, Delft University of Technology, 2004-2013.
