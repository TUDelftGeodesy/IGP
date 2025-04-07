function savepng(h,pltdir)
%SAVEPNG   Save figure(s) as png file(s).
%  SAVEPNG(HP,PLTDIR) saves figures with plot handles in HP
%  to a png file in the directory PLTDIR. 

% Created:   3 April 2019 by Hans van der Marel
% Modified:

if nargin < 2
    pltdir='';
elseif ~isdir(pltdir)
    mkdir(pltdir);
end

for k=1:numel(h)
   set(h(k),'PaperPositionMode','Auto');
   fprintf('saving %s.png\n',fullfile(pltdir,h(k).Name))
   print(h(k),'-dpng',fullfile(pltdir,h(k).Name));
end

end