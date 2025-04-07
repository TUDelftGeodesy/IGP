function h=pltarrow(from,to,varargin)
% pltarrow   plot an arrow
%    PLTARROW(FROM,TO) plots an arrow from position FROM to position TO
%
%    PLTARROW(FROM,TO,...) other arguments will be passed to the plot function.
%
%    H=PLTARROW(...) returns the plot handle H.
%
%    See also plot. 

if nargin < 2 
   error('not the correct number of input arguments')
end
if isempty(varargin)
  varargin{1}='b';
end

% Arrow head parameters

alpha = 0.33;   % Size of arrow head relative to the length of the vector
beta = 0.33;    % Width of the base of the arrow head relative to the length

UpArrow= [ -alpha 0 alpha  ; ...      
           -beta  0 -beta  ];
scale=1;

% Get coordinates

x1=from(1); y1=from(2);
x2=to(1); y2=to(2);

% Scale and rotate the arrow head

u=x2-x1;
v=y2-y1;
th=-pi/2+atan2(v,u);
R=[cos(th) -sin(th); sin(th) cos(th)];   % Rotation matrix for local slope of vv.
A=scale*R*UpArrow;    % Rotate the arrowhead.

% Plot the arrow

hold_state = get(gca,'nextplot');
hold on
patch(x2+A(1,:),y2+A(2,:),varargin{:});
plot(x2+A(1,:),y2+A(2,:),varargin{:});
h=plot([x1 x2],[y1 y2],varargin{:});
get(h,'MarkerFaceColor');
get(h,'MarkerEdgeColor');
set(gca,'nextplot',hold_state);

end