function hf=tseriesplot(tseries,varargin)
%TSERIESPLOT  Plot time series components.
%   TSERIESPLOT(TSERIES) plots timeseries fitted components from the time 
%   series structure TSERIES.
%
%   TSERIESPLOT(TSERIES,OPTNAME,OPTVALUE,...) specify option name value pairs
%      doplots              plot options: 
%                scalar        0:none, 1: basic, 2:extensive, if negative, visible is off
%                cell array    {'fit' 'comp' 'resid' 'stdev' 'neu' 'dneu' }
%      saveplot             if true save the plots to disk
%
%   HF=TSERIESPLOT(...) returns the plot handles.
%
%   See also tseriesfit.
%
%   (c) Hans van der Marel, Delft University of Technology, 2014-2016.

% Default options

opt.doplots=0;               % by default we don't do plotting, 1 is basic plotting, 2 extensive
opt.visible='on';
opt.saveplots=false;         % by default plots are not saved 
opt.plotdir='';              % directory for plots 
opt.empstdev=true;           % scale st.dev. with emperical factor
opt.id='';                   % optional id for title (will be overwritten by series id)

% Check the input arguments

if ~isstruct(tseries)
   error('First input argument must be a tseries structure.')
end
for k=1:2:numel(varargin)-1
   if isfield(opt,varargin{k})
      opt.(varargin{k})=varargin{k+1};
   else
      error(['Illegal option ' varargin{k} ])
   end
end

% Check legacy plot options

if isnumeric(opt.doplots)
  if opt.doplots < 0
     visible='off';
     opt.doplots=abs(opt.doplots);
  else
     visible='on';
  end
  if opt.doplots > 1
    doplots={'fit' 'comp' 'resid' 'stdev' 'neu' 'dneu' };
  elseif opt.doplots > 0
    doplots={'fit' 'comp' 'resid' }; 
  else
    hf=[];
    return;
  end
elseif ischar(opt.doplots)
  doplots={opt.doplots};
  visible=opt.visible;    
else
  doplots=opt.doplots;
  visible=opt.visible;
end

% Get the station name, epoch, neu and sneu data

station=tseries.station;
if isfield(tseries,'year')
  year=tseries.year;
else
  epoch=tseries.epoch;
  year=date2dyear(epoch);
end
neu=tseries.neu;
sneu=tseries.sneu;

% Get solution id

if isfield(tseries,'id')
  id=tseries.id;
else
  id=opt.id;
end
savename='';
if ~isempty(id)
  id=[' ' id];
  savename=['_' id];
end

% Get antenna and receiver names and other events

%  EVENTS is a structure array with at least the following fields
%    events(k).year   start of the k'th event in decimal years
%    events(k).type   type of event: REC=receiver change, ANT=antenna
%                     change, or others.
%    events(k).name   name of the receiver or antenna (only for REC or ANT)

if isfield(tseries,'events')
   events=tseries.events;
else
   events=[];
end

% Get fitted components

neuresidual=tseries.neuresidual;
neunojumps=tseries.neu-tseries.neujumps;
neufit=tseries.neufit-tseries.neujumps;
neutrend=tseries.neutrend;
neutrendldn=tseries.neutrend+tseries.neuatmld;

% rename dneu and neu (for historic purposes

dneu=neu;
neu=dneu+(year-tseries.t0)*tseries.vneu0;    

% Compute scaling factors for standard deviations

empstd=std(diff(dneu))./sqrt(2)*1000;
meanstd=mean(sneu)*1000;

f=empstd./meanstd;

fprintf('\n                    North [mm]  East [mm]    Up [mm]\n')
fprintf('Emperical St.Dev.:  %10.3f %10.3f %10.3f\n',empstd);
fprintf('Formal St.Dev.:     %10.3f %10.3f %10.3f\n',meanstd);
fprintf('Factor (Estimated): %10.3f %10.3f %10.3f\n',f);
if opt.empstdev
  fprintf('\nStandard deviations will be re-scaled.\n') 
else
  fprintf('\nStandard deviations will NOT be re-scaled.\n') 
  f=[ 1 1 1];
end  

% Get screensize and position figure [ left bottom width height ]

scrsz = get(0,'ScreenSize');
%figpos=[ scrsz(1)+scrsz(3)*0.05 scrsz(2)+scrsz(4)*0.05 scrsz(3)*0.9 scrsz(4)*0.9]; 
figpos=[ scrsz(3)*0.05 scrsz(4)*0.05 scrsz(3)*0.9 scrsz(4)*0.94]; 

% Do the plotting

for k=1:numel(doplots)
  selplot=lower(doplots{k});
  switch selplot

     case 'neu'
         
        hf(k)=figure('Name',[ station id ' (IGS08b)' ] ,'NumberTitle','off','OuterPosition',figpos,'visible',visible);

        subplot(3,1,1)
        pltseries1b(year,neu(:,1),sneu(:,1)*f(1),events)
        ylabel('\Delta N [mm]')
        title([ station id ' (IGS08b)' ]  )
        subplot(3,1,2)
        pltseries1b(year,neu(:,2),sneu(:,2)*f(2),events)
        ylabel('\Delta E [mm]')
        subplot(3,1,3)
        pltseries1b(year,neu(:,3),sneu(:,3)*f(3),events)
        ylabel('\Delta U [mm]')

        if opt.saveplots
           print('-dpng','-r 300',fullfile(opt.plotdir,[ station savename '_IGS08b.png' ]))
           %print('-depsc',[ 'plots/' station savename '_IGS08b.eps' ])
        end

     case 'dneu' 

        hf(k)=figure('Name',[ station id ' (ETRS89)' ] ,'NumberTitle','off','OuterPosition',figpos,'visible',visible);

        subplot(3,1,1)
        pltseries1b(year,dneu(:,1),sneu(:,1)*f(1),events,[-12 12])
        ylabel('\Delta N [mm]')
        title([ station id ' (ETRS89)' ]  )
        subplot(3,1,2)
        pltseries1b(year,dneu(:,2),sneu(:,2)*f(2),events,[-12 12])
        ylabel('\Delta E [mm]')
        subplot(3,1,3)
        pltseries1b(year,dneu(:,3),sneu(:,3)*f(3),events)
        ylabel('\Delta U [mm]')

        if opt.saveplots
           print('-dpng','-r 300',fullfile(opt.plotdir,[ station savename '_ETRS89.png' ]))
           %print('-depsc',[ 'plots/' station savename '_ETRS89.eps' ])
        end

     case {'series','fit'} 

        hf(k)=figure('Name',[ station id ' (Series w/ fit)' ] ,'NumberTitle','off','OuterPosition',figpos,'visible',visible);

        subplot(3,1,1)
        pltseries1b(year,neunojumps(:,1),sneu(:,1)*f(1),events,[-12 12])
        h = findobj(gca,'Type','patch');
        h1 = findobj(gca,'Type','line','-and','Color','blue');
        h2=plot(year,neufit(:,1)*1000,'g-');    
        h3=plot(year,neutrendldn(:,1)*1000,'r-');    
        h4=plot(year,neutrend(:,1)*1000,'k-');    
        ylabel('\Delta N [mm]')
        title([ station id ' (Series w/ Fit)' ]  )
        hl=legend([ h1(1) h2 h4 h3 h(1)],'Obs','Fit','Trend','Trend + Atm.Ld.','2 sigma','Location','SouthWest','Orientation','horizontal');
        plabel=get(hl,'Position');
        plabel(1)=0.5-plabel(3)/2;
        plabel(2)=0.02;
        set(hl,'Position',plabel);
        subplot(3,1,2)
        pltseries1b(year,neunojumps(:,2),sneu(:,2)*f(2),events,[-12 12])
        plot(year,neufit(:,2)*1000,'g-')    
        plot(year,neutrendldn(:,2)*1000,'r-')    
        plot(year,neutrend(:,2)*1000,'k-')    
        ylabel('\Delta E [mm]')
        subplot(3,1,3)
        pltseries1b(year,neunojumps(:,3),sneu(:,3)*f(3),events,[-12 12])
        plot(year,neufit(:,3)*1000,'g-')    
        plot(year,neutrendldn(:,3)*1000,'r-')    
        plot(year,neutrend(:,3)*1000,'k-')    
        ylabel('\Delta U [mm]')

        if opt.saveplots
           print('-dpng','-r 300',fullfile(opt.plotdir,[ station savename '_series.png' ]))
           %print('-depsc',[ 'plots/' station savename '_series.eps' ])
        end

     case 'comp' 

        hf(k)=figure('Name',[ station id ' (Signal components)' ] ,'NumberTitle','off','OuterPosition',figpos,'visible',visible);

        subplot(3,1,1)
        pltseries1b(year,[],sneu(:,1)*f(1),events,[-6 6]);
        h = findobj(gca,'Type','patch');
        h2=plot(tseries.year,tseries.neutempi(:,1)*1000,'g-');    
        h23=plot(tseries.year,tseries.neuharmonic(:,1)*1000+tseries.neutempi(:,1)*1000,'b-');    
        h3=plot(tseries.year,tseries.neuharmonic(:,1)*1000,'k-');    
        h4=plot(tseries.year,tseries.neuatmld(:,1)*1000,'r-');    
        ylabel('\Delta N [mm]')
        title([ station id ' (Periodic + Temp.Influence + Atm.Loading)' ]  )
        hl=legend([ h2 h3 h23 h4 h(1)],'Temp.Infl.','Periodic','Temp+Periodic','Atm.Ld.','2 sigma','Location','SouthWest','Orientation','horizontal');
        plabel=get(hl,'Position');
        plabel(1)=0.5-plabel(3)/2;
        plabel(2)=0.02;
        set(hl,'Position',plabel);
        subplot(3,1,2)
        pltseries1b(year,[],sneu(:,2)*f(2),events,[-6 6])
        plot(tseries.year,tseries.neutempi(:,2)*1000,'g-')    
        plot(tseries.year,tseries.neuharmonic(:,2)*1000+tseries.neutempi(:,2)*1000,'b-');    
        plot(tseries.year,tseries.neuharmonic(:,2)*1000,'k-')    
        plot(tseries.year,tseries.neuatmld(:,2)*1000,'r-')    
        ylabel('\Delta E [mm]')
        subplot(3,1,3)
        pltseries1b(year,[],sneu(:,3)*f(3),events,[-6 6])
        plot(tseries.year,tseries.neutempi(:,3)*1000,'g-')    
        plot(tseries.year,tseries.neuharmonic(:,3)*1000+tseries.neutempi(:,3)*1000,'b-');    
        plot(tseries.year,tseries.neuharmonic(:,3)*1000,'k-')    
        plot(tseries.year,tseries.neuatmld(:,3)*1000,'r-')    
        ylabel('\Delta U [mm]')

        if opt.saveplots
           print('-dpng','-r 300',fullfile(opt.plotdir,[ station savename '_components.png' ]))
           %print('-depsc',[ 'plots/' station savename '_components.eps' ])
        end
        
     case 'resid'

        hf(k)=figure('Name',[ station id ' (Residuals)' ] ,'NumberTitle','off','OuterPosition',figpos,'visible',visible);

        subplot(3,1,1)
        pltseries1b(year,neuresidual(:,1),sneu(:,1)*f(1),events,[-12 12])
        ylabel('\Delta N [mm]')
        title([ station id ' (Residuals)' ]  )
        subplot(3,1,2)
        pltseries1b(year,neuresidual(:,2),sneu(:,2)*f(2),events,[-12 12])
        ylabel('\Delta E [mm]')
        subplot(3,1,3)
        pltseries1b(year,neuresidual(:,3),sneu(:,3)*f(3),events,[-12 12])
        ylabel('\Delta U [mm]')

        if opt.saveplots
           print('-dpng','-r 300',fullfile(opt.plotdir,[ station savename '_residuals.png' ]))
           %print('-depsc',[ 'plots/' station savename '_residuals.eps' ])
        end


     case 'stdev'

        hf(k)=figure('Name',[ station id ' (St.Dev.)' ] ,'NumberTitle','off','OuterPosition',figpos,'visible',visible);

        subplot(3,1,1)
        pltseries1b(year,sneu(:,1)*f(1),[],events)
        ylabel('\sigma N [mm]')
        title([ station id ' (St.Dev.)' ]  )
        subplot(3,1,2)
        pltseries1b(year,sneu(:,2)*f(2),[],events)
        ylabel('\sigma E [mm]')
        subplot(3,1,3)
        pltseries1b(year,sneu(:,3)*f(3),[],events)
        ylabel('\sigma U [mm]')

        if opt.saveplots
           print('-dpng','-r 300',fullfile(opt.plotdir,[ station savename '_stdev.png' ]))
           %print('-depsc',[ 'plots/' station savename '_stdev.eps' ])
        end

      otherwise
          
        printf('Unknown plot type %s\n',selplot);         

  end
  
end

end

