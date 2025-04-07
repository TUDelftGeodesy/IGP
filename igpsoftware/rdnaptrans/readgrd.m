function gstruct=readgrd(filename)
%READGRD  Read binary grid files in the format of the program Surfer(R) 
%  Read binary grid files in the format of the program Surfer(R) 
%  Input:
%    filename      name of the binary file
%  Output:
%    gstruct       structure with grid data
%      gstruct.x     vector with the x coordinates (longitude)
%      gstruct.y     vector with the y coordinates (latitude)
%      gstruct.mat   matrix with values of dimension [length(x),length(y)]
%
%  Layout of mat:
%
%    +--> y                                        West 
%    |                                              ^
%    v   (1,1)  (1,2) .. (1,ny)                     |
%        (2,1)  (2,2)                       South <-+-> North
%    x     :                                        |
%        (nx,1)          (nx,ny)                    v
%                                                  East
%
%  See also grdint and testgrdint.
%
%  (c) Hans van der Marel, Delft University of Technology, 2004.

% Created:   6 Dec 2004 by Hans van der Marel, TUD
% Modified: 

% Open grid file

fid = fopen(filename,'r');

if (fid < 0) 
  disp(['Error opening file ',filename]);
  return;
end 

% Read file id

id = char(fread(fid,4))';

if (id ~= 'DSBB') 
  disp([filename,' is not a valid grd file']);
  fclose(fid);
  return;
end 

%  Read header parameters

nx=fread(fid,1,'int16');
ny=fread(fid,1,'int16');

minx=fread(fid,1,'double');
maxx=fread(fid,1,'double');
miny=fread(fid,1,'double');
maxy=fread(fid,1,'double');
minz=fread(fid,1,'double');
maxz=fread(fid,1,'double');

dx = (maxx-minx)/(nx-1);
dy = (maxy-miny)/(ny-1);

gstruct.x=minx:dx:maxx;
gstruct.y=miny:dy:maxy;

% Read body of the file

gstruct.mat=fread(fid,[nx,ny],'float32');

% Replace illegal values by NaN

gstruct.mat(gstruct.mat > maxz | gstruct.mat < minz)=NaN;

% Done, close file

fclose(fid);

% Some plots for debugging

%figure; contourf(x,y,mat');axis equal; axis tight;colorbar;
%figure; imagesc(x,y,mat');axis equal;axis tight;set(gca,'Ydir','Normal');colorbar;

return
