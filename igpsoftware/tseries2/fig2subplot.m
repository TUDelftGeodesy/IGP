function ht=fig2subplot(hp)
%FIG2SUBPLOT   Copy figures to subplots.
%  HT=FIG2SUBPLOT(HP) copies figures with plot handles in the matrix HP
%  to a new figure of subplots with handle HT. The shape of matrix HP
%  with figure handles defines the subplot layout.

% Created:  3 February 2019 by Hans van der Marel
% Modified:

% Set subplot margings

mx=[.03 .01];
my=[.04 .04];

% Create a (almost) full screen figure

ht=figure('units','normalized','outerposition',[0 .05 1 .95]);

% Copy the figures with plot handles in HP to destination figure

[m,n]=size(hp);
for i=1:m
  for j=1:n
    if isempty(hp(i,j)), continue;, end
    % Now copy contents of each figure over to destination figure
    % Modify position of each axes as it is transferred
    try
      h = get(hp(i,j),'Children');
      newh = copyobj(h,ht);
      for k = 1:length(newh)
         set(newh(k),'Position',...
         [ (j-1)/n+mx(1) (m-i)/m+my(1) 1/n-mx(1)-mx(2) 1/m-my(1)-my(2)]);
      end
    end
  end
end

end