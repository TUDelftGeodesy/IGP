function mstruct=nlbessel()
%NLBESSEL Parameters for the Dutch Stereographic double projection 
%  Parameters for the Dutch Stereographic double projection 
%  Syntax
%            mstruct=nlbessel()
%  Input
%    none
%  Output
%    mstruct Structure with projection parameters
%              Geographic NL-Bessel coordinates of Amersfoort (pivot point and 
%              projection base point)
%                mstruct.origin(1)  latitude in radians
%                mstruct.origin(2)  longitude in radians
%              Parameters of ellipsoids Bessel1841 
%                mstruct.a          half major axis in meters
%                mstruct.inv_f      inverse flattening
%              Parameters of RD projection
%                mstruct.scale      scale factor (k in some notations)
%                                   this factor was first defined by Hk.J. Heuvelink as 
%                                   pow(10,-400.00e-7), nowadays we define it as exactly 0.9999079
%                mstruct.falseeasting
%                mstruct.falsenorthing 
%
%  References:
%    G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning 
%      at Sea". Delft University of Technology, 1995.
%    G. Strang van Hees, "Globale en lokale geodetische systemen", 
%      Nederlandse Commissie voor Geodesie (NCG), Delft, 1997.
%    Hk.J. Heuvelink, "De stereografische kaartprojectie in hare toepassing bij de Rijksdriehoeksmeting". Delft: Rijkscommissie voor Graadmeting en Waterpassing, 1918. 
%    HTW, "Handleiding voor de Technische Werkzaamheden van het Kadaster". Apeldoorn: Kadaster, 1996.
%
%  See also stereodp.
%
%  (c) Hans van der Marel, Delft University of Technology, 2004.

% Created:   6 Dec 2004 by Hans van der Marel, TUD
% Modified:

PHI_AMERSFOORT_BESSEL    = ( 52.0+ 9.0/60.0+22.178/3600.0 ) ;
LAMBDA_AMERSFOORT_BESSEL = (  5.0+23.0/60.0+15.500/3600.0 ) ;
H_AMERSFOORT_BESSEL      =  0.0;

A_BESSEL     = 6377397.155;   
INV_F_BESSEL = 299.1528128;

SCALE_RD = 0.9999079;
X_AMERSFOORT_RD = 155000;
Y_AMERSFOORT_RD = 463000;

mstruct.origin(1)=PHI_AMERSFOORT_BESSEL;
mstruct.origin(2)=LAMBDA_AMERSFOORT_BESSEL;
mstruct.a=A_BESSEL;
mstruct.inv_f=INV_F_BESSEL;
mstruct.scale=SCALE_RD;
mstruct.falseeasting=X_AMERSFOORT_RD;
mstruct.falsenorthing=Y_AMERSFOORT_RD;

return;
