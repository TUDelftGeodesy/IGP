
%%
pntcrd = rand([2000,2])*100;

%% BBOX or Polygon [OK}

bbox = [20 20 ; 50 50 ]
poly = roi2poly(bbox)

pntmask = getpntmask(pntcrd,bbox);
plotmask(pntcrd,pntmask)

pntmask = getpntmask(pntcrd,poly);
plotmask(pntcrd,pntmask)


%% Two Plygons [OK]

poly1 = roi2poly([20 20 ; 50 50 ]);
poly2 = roi2poly([70 20 ; 90 90 ]);

poly = [ poly1 ; NaN NaN ; poly2 ];
pntmask = getpntmask(pntcrd,poly);

plotmask(pntcrd,pntmask)

%% Two Polygons - cw for hole [OK]

poly1 = roi2poly([20 20 ; 60 70 ]);
poly2 = roi2poly([30 30 ; 40 50 ]);

[x, y] = poly2cw(poly2(:,1),poly2(:,2));
poly2 = [ x y ];

poly = [ poly1 ; NaN NaN ; poly2 ];
pntmask = getpntmask(pntcrd,poly);

plotmask(pntcrd,pntmask)


%% Shape file with two polygons - two ways [OK]

poly1 = roi2poly([20 20 ; 50 50 ]);
poly2 = roi2poly([70 20 ; 90 90 ]);

S = geoshape(poly1(:,1),poly1(:,2),'Geometry','Polygon');
S(2)=poly2;

S

shapewrite(S, 'test_test.shp')

poly=roi2poly('test_test.shp')

pntmask = getpntmask(pntcrd,'test_test.shp');
plotmask(pntcrd,pntmask)

poly = [ poly1 ; NaN NaN ; poly2 ];

S = geoshape(poly(:,1),poly(:,2),'Geometry','Polygon');

S

shapewrite(S, 'test_test.shp')

poly=roi2poly('test_test.shp')

pntmask = getpntmask(pntcrd,'test_test.shp');
plotmask(pntcrd,pntmask)


%% KML file with two polygons - two ways [OK]

poly1 = roi2poly([20 20 ; 50 50 ]);
poly2 = roi2poly([70 20 ; 90 90 ]);

S = geoshape(poly1(:,1),poly1(:,2),'Geometry','Polygon');
S(2)=poly2;

S

kmlwrite('test_test1.kml',S)

poly=roi2poly('test_test1.kml')

pntmask = getpntmask(pntcrd,'test_test1.kml');
plotmask(pntcrd,pntmask)

poly = [ poly1 ; NaN NaN ; poly2 ];

S = geoshape(poly(:,1),poly(:,2),'Geometry','Polygon');

S

kmlwrite('test_test2.kml',S)

poly=roi2poly('test_test2.kml')

pntmask = getpntmask(pntcrd,'test_test2.kml');
plotmask(pntcrd,pntmask)

%% shapefile with holes - two ways [OK]

poly1 = roi2poly([20 20 ; 60 70 ]);
poly2 = roi2poly([30 30 ; 40 50 ]);

[x, y] = poly2cw(poly2(:,1),poly2(:,2));
poly2 = [ x y ];

S = geoshape(poly1(:,1),poly1(:,2),'Geometry','Polygon');
S(2)=poly2;

S

shapewrite(S, 'test_test.shp')

poly=roi2poly('test_test.shp')

pntmask = getpntmask(pntcrd,'test_test.shp');
plotmask(pntcrd,pntmask)

poly = [ poly1 ; NaN NaN ; poly2 ];

S = geoshape(poly(:,1),poly(:,2),'Geometry','Polygon');

S

shapewrite(S, 'test_test.shp')

poly=roi2poly('test_test.shp')


pntmask = getpntmask(pntcrd,'test_test.shp');
plotmask(pntcrd,pntmask)

%% KML with holes - two ways [FAILS]

poly1 = roi2poly([20 20 ; 60 70 ]);
poly2 = roi2poly([30 30 ; 40 50 ]);

[x, y] = poly2cw(poly2(:,1),poly2(:,2));
poly2 = [ x y ];

S = geoshape(poly1(:,1),poly1(:,2),'Geometry','Polygon');
S(2)=poly2;

S

kmlwrite('test_test3.kml',S)

poly=roi2poly('test_test3.kml')

pntmask = getpntmask(pntcrd,'test_test3.kml');
plotmask(pntcrd,pntmask)

poly = [ poly1 ; NaN NaN ; poly2 ];

S = geoshape(poly(:,1),poly(:,2),'Geometry','Polygon');

S

kmlwrite('test_test4.kml',S)

poly=roi2poly('test_test4.kml')


pntmask = getpntmask(pntcrd,'test_test4.kml');
plotmask(pntcrd,pntmask)

%%

function plotmask(pntcrd,pntmask)
figure
plot(pntcrd(pntmask,1),pntcrd(pntmask,2),'r*')
hold on
plot(pntcrd(~pntmask,1),pntcrd(~pntmask,2),'b.')
end

