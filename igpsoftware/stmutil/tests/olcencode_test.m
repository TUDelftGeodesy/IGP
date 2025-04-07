
% For the last two elements we expect '9F4649CQ+HC' and '6PH57VP3+PR6' (one
% lsd added) 7FG49QCJ+2VXGJ 

lat = [ 52.3 45.6 23.4254235 -45.0 0.  52.121437   1.286785    20.3701135 52.121437  20.3701135 52.121437]
lon = [  4    6   -130.4212  130.5222  12 4.388562 103.854503 2.78223535156 4.388562  2.78223535156 4.388562]

olc=olcencode(lat,lon)

olc=olcencode(lat,lon,5)
olc=olcencode(lat,lon,1)

olc=olcencode(lat(1:3),lon(1:3))

olc=olcencode(lat(1),lon(1))