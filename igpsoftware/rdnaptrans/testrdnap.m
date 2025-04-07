function testrdnap(varargin)
%TESTRDNAP  Test RDNAPTRANS using a series of test points.
%  This funnction tests the RDNNAPTRANS transformation using a series of test
%  points with given coordinates in ETRS89 and RDNAP.
%  Syntax
%             testrdnap()
%             testrdnap(options)
%  Options
%    id       system [RDNAPTRANS2008|RDNAPTRANS2004], optional, default most 
%             recent
%    extrd    RD correction grid extension [NaN|zero], default NaN
%  Output
%    none
%
%  See also rdnap2etrs and etrs2rdnap.
% 
%  (c) Hans van der Marel, Delft University of Technology, 2004-2013

% Created:   6 Dec 2004 by Hans van der Marel, TUD
% Modified:  6 Jun 2013 by Hans van der Marel
%              - upgraded to a function
%              - Select between RDNAPTRANS2004 and RDNAPTRANS2008 
%              - Adapted to use new 'plh' option'

% if scalar==1 test the scalar calls instead of vectorized call

scalar=0;

% process the options

id='MOSTRECENT';
extrd='NAN';

for i=1:length(varargin)
  switch upper(varargin{i})
     case {'RDNAPTRANS2008','2008','RDNAPTRANS2004','2004','MOSTRECENT','CURRENT'}
       id=varargin{i};
     case {'NAN','ZERO'}
       extrd=varargin{i};
     otherwise
       error('invalid option')
  end
end

% Test coordinates for RDNAPTRANS(tm)2004/2008 
% --------------------------------------------------------------------------
% The RDNAPTRANS(tm)2004/2008 procedure should be tested in both directions. The 
% differences with the coordinates listed below should be smaller or equal to:
% RD x and y coordinates:                      0.0010 meters
% NAP heights and ETRS89 ellipsoidal heights:  0.0010 meters
% ETRS89 latitude and longitude:               0.0000000100 degrees

names=[
'Texel             ';
'Noord-Groningen   ';
'Amersfoort        ';
'Amersfoort_100m   ';
'Zeeuws-Vlaanderen ';
'Zuid-Limburg      ';
'Maasvlakte        ';
'outside           ';
'no_rd_geoid       ';
'no_geoid          ';
'no_rd             ';
'edge_rd           '];

rdnap=[
   117380.1200 575040.3400   1.0000  ;
   247380.5600 604580.7800   2.0000  ;
   155000.0000 463000.0000   0.0000  ;
   155000.0000 463000.0000 100.0000  ;
    16460.9100 377380.2300   3.0000  ;
   182260.4500 311480.6700 200.0000  ;
    64640.8900 440700.0100   4.0000  ;
   400000.2300 100000.4500   5.0000  ;
   100000.6700 300000.8900   6.0000  ;
   100000.6700 350000.8900   6.0000  ;
    79000.0100 500000.2300   7.0000  ;
    50000.4500 335999.6700   8.0000  ];

switch upper(id)
   case {'RDNAPTRANS2008','2008','MOSTRECENT','CURRENT'}
     plhetrs=[
       53.160753042      4.824761912     42.8614     ;
       53.419482050      6.776726674     42.3586     ;
       52.155172897      5.387203657     43.2551     ;
       52.155172910      5.387203658     143.2551    ;
       51.368607152      3.397588595     47.4024     ;
       50.792584916      5.773795548    245.9478     ;
       51.947393898      4.072887101     47.5968     ;
       48.843030210      8.723260235     52.0289     ;
       50.687420392      4.608971813     51.6108     ;
       51.136825197      4.601375361     50.9672     ;
       52.482440839      4.268403889     49.9436     ;
       51.003976532      3.891247830     52.7427     ];
   case {'RDNAPTRANS2004','2004'}
     plhetrs=[
       53.1607530501     4.8247619722      42.8702   ;
       53.4194820671     6.7767267456      42.3674   ;
       52.1551728997     5.3872037291      43.2639   ;
       52.1551729135     5.3872037298     143.2639   ;
       51.3686071451     3.3975886580      47.4112   ;
       50.7925849091     5.7737956320      245.9566  ;
       51.9473938958     4.0728871651      47.6056   ;
       49.0000000000     9.0000000000      43.0000   ;
       50.6874122400     4.6089625000      44.0000   ;
       51.1368171300     4.6013663100      45.0000   ;
       52.4824386900     4.2684038100      46.0000   ;
       51.0039764500     3.8912478600      47.0000   ];
   otherwise
     error(['Unknown system ' id ])
end

plhetrs=[plhetrs(:,1:2)*pi/180 plhetrs(:,3)];

% From ETRS89 to RD/NAP

fprintf('ETRS89 to RD/NAP discrepancies:\n\n');
fprintf('name                     dx (m)       dy (m)       dH (m)\n\n');
if ~scalar
   rdnap1=etrs2rdnap(plhetrs,varargin{:},'plh');
end
for i=1:length(plhetrs)
   if scalar
      rdnap1(i,1:3)=etrs2rdnap(plhetrs(i,1:3),varargin{:},'plh');
   end
   fprintf('%s %12.6f %12.6f %12.6f\n',names(i,:),rdnap1(i,:)-rdnap(i,:))
end
fprintf('\n\n');

% From RD/NAP to ETRS

fprintf('Inverse transformation RD/NAP->ETRS:\n\n');
fprintf('name                 dPhi (mas)   dLam (mas)       dh (m)\n\n');
if ~scalar
   plhetrs1=rdnap2etrs(rdnap1,varargin{:},'plh');
end
for i=1:length(rdnap1)
   %if any(isnan(rdnap1(i,:))), break;, end
   if scalar
      plhetrs1(i,1:3)=rdnap2etrs(rdnap1(i,1:3),varargin{:},'plh');
   end
   detrs=plhetrs1(i,:)-plhetrs(i,:);
   detrs(1:2)=detrs(1:2)*180/pi*60*60*1000;
   fprintf('%s %12.6f %12.6f %12.6f\n',names(i,:),detrs)
end
fprintf('\n\n');

% From RD/NAP to ETRS

fprintf('RD/NAP to ETRS89 discrepancies:\n\n');
fprintf('name                 dPhi (mas)   dLam (mas)       dh (m)\n\n');
if ~scalar
   plhetrs2=rdnap2etrs(rdnap,varargin{:},'plh');
end
for i=1:length(rdnap)
   if scalar
      plhetrs2(i,1:3)=rdnap2etrs(rdnap(i,1:3),varargin{:},'plh');
   end
   detrs=plhetrs2(i,:)-plhetrs(i,:);
   detrs(1:2)=detrs(1:2)*180/pi*60*60*1000;
   fprintf('%s %12.6f %12.6f %12.6f\n',names(i,:),detrs)
end
fprintf('\n\n');

% From RD/NAP to ETRS inverse transformation

fprintf('Inverse transformation ETRS->RD/NAP:\n\n');
fprintf('name                     dx (m)       dy (m)       dH (m)\n\n');
if ~scalar
   rdnap2=etrs2rdnap(plhetrs2,varargin{:},'plh');
end
for i=1:length(plhetrs2)
   %if any(isnan(plhetrs2(i,:))), break;, end
   if scalar
      rdnap2(i,1:3)=etrs2rdnap(plhetrs2(i,1:3),varargin{:},'plh');
   end
   fprintf('%s %12.6f %12.6f %12.6f\n',names(i,:),rdnap2(i,:)-rdnap(i,:))
end
fprintf('\n\n');

return