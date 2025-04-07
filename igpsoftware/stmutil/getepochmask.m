function epochMask=getepochmask(epochDyears,poi) 
%GETEPOCHMASK  Compute epoch mask for Period Of Interest (POI). 
%   EPOCHMASK=GETEPOCHMASK(EPOCHS,POI) compute the epoch mask EPOCHMASK 
%   for the array EPOCHS with timestamps (dyear, matlab datenumber, ...)  
%   for the Period of Interest POI. The Period of Interest POI is a
%   Matlab array with two columns, with respectively start and end
%   times. The function does not check for overlaps.
%
%   Examples:
%
%      epochMask = getepochmask(epochDyear,[ 2010 2021.5])
%      epochMask = getepochmask(epochDyear,[ 2010 2012.2 ; 2015.3 2021.5])
%
%   (c) Hans van der Marel, Delft University of Technology, 2021. 

% Created:  22 Oct 2021 by Hans van der Marel
% Modified: 
%

if isempty(poi)
   epochMask = true(size(epochDyears));
else
   epochMask = false(size(epochDyears));
   for k=1:size(poi,1)
      epochMask = epochMask | ( epochDyears >= poi(k,1) & epochDyears <= poi(k,2) );
   end
end

end