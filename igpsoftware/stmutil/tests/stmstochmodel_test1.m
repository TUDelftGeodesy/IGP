%% Test script for stmstochmodel - Basic tests
%
% Hans van der Marel
%
%% Simulate coordinates [km] and time [dyears]

%[x,y]=meshgrid(-2:.1:2,-2:.1:2);
%xy=[ x(:) y(:)];

numPoints=7;   
xy=[ -2+rand(numPoints,2)*4 ];

t=2001:2010;
t=t+rand(size(t))*.5+.3;

%% Test temporal models 

C=stmstochmodel('PL(alpha=0,sigma=0.5,fs=1,ndays=1)',[],[],t)
C=stmstochmodel('PL(alpha=-0.8,sigma=0.5,fs=1,ndays=1)',[],[],t)
C=stmstochmodel('PL(alpha=-1,sigma=0.5,fs=1,ndays=1)',[],[],t)
C=stmstochmodel('PL(alpha=-2,sigma=0.5,fs=1,ndays=1)',[],[],t)

%%
C=stmstochmodel('PL(alpha=-1,sigma=0.5,fs=1,ndays=1)',[],[],t)
C=stmstochmodel('PL(alpha=-1,sigma=0.5,fs=1,ndays=2)',[],[],t)
C=stmstochmodel('PL(alpha=-1,sigma=0.5,fs=1,ndays=3)',[],[],t)
C=stmstochmodel('PL(alpha=-1,sigma=0.5,fs=1,ndays=10)',[],[],t)

%%
C=stmstochmodel('iGGM(alphaR=-1,alphaL=-.8,gamma=0.99,sigma=0.5,fs=1,ndays=1)',[],[],t)
C=stmstochmodel('iGGM(alphaR=-1,alphaL=-.8,gamma=0.99,sigma=0.5,fs=1,ndays=2)',[],[],t)
C=stmstochmodel('iGGM(alphaR=-1,alphaL=-.8,gamma=0.99,sigma=0.5,fs=1,ndays=5)',[],[],t)

%%
C=stmstochmodel('FOGM(rho=1/365,sigma=1.5,fs=1,ndays=1)',[],[],t)
figure;imagesc(C);colorbar;
C=stmstochmodel('FOGM(rho=1/365,sigma=1.5,fs=1,ndays=4)',[],[],t)
figure;imagesc(C);colorbar;

%%
epochAttrib.ndays=ones(10,1)*2
%C=stmstochmodel('PL(alpha=-1,sigma=0.5,fs=1,ndays=ndays)',[],[],t)
%C=stmstochmodel('PL(alpha=-1,sigma=0.5,fs=1,ndays=wrong)',[],[],t,[],epochAttrib)
C=stmstochmodel('PL(alpha=-1,sigma=0.5,fs=1,ndays=ndays)',[],[],t,[],epochAttrib)

%%
C=stmstochmodel({'iGGM(alphaR=-1,alphaL=-.8,gamma=0.99,sigma=0.5,fs=1,ndays=1)','PL(alpha=-1,sigma=0.5,fs=1,ndays=2)'},[],[],t);
figure;imagesc(C);colorbar;

%% Test temporal models with spatial correlation

C=stmstochmodel({ ...
    'iGGM(alphaR=-1,alphaL=-.8,gamma=0.99,sigma=0.5,fs=1,ndays=1)', ...
    'PL(alpha=-1,sigma=0.5,fs=1,ndays=2)', ...
    'spatialcov(rho=.88)' ...
    },[],xy,t);
figure;imagesc(C);colorbar;

C=stmstochmodel({ ...
    'iGGM(alphaR=-1,alphaL=-.8,gamma=0.99,sigma=0.5,fs=1,ndays=1)', ...
    'PL(alpha=-1,sigma=0.5,fs=1,ndays=2)', ...
    'spatialcov(rho=.88)', ...
    'WN(sigma=3)' ...
    },[],xy,t);
figure;imagesc(C);colorbar;


%% Test white noise models without spatial correlation

C=stmstochmodel({ ...
    'WN(sigma=2)' ...
    'spatialcov(rho=0)' ...
    },[],xy,t);
figure;imagesc(C);colorbar;

C=stmstochmodel({ ...
    'WN(sigma=2)' ...
    'spatialcov(rho=Inf)' ...
    },[],xy,t);
figure;imagesc(C);colorbar;

C=stmstochmodel({ ...
    'WN(sigma=2)' ...
    },[],xy,t);
figure;imagesc(C);colorbar;

C=stmstochmodel({ ...
    'WN(sigma=2)' ...
    },[],[],t);
figure;imagesc(C);colorbar;



%% Test InSAR models

C=stmstochmodel('tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)', ...
          [],xy,t);                                          % Radarsat-2
figure;imagesc(C);colorbar;

%%      

mask=true(size(xy,1),numel(t));
mask(5:6,3:6)=false;
C=stmstochmodel('tudinsar4(s20=9.49,s2t=4.53,s2s=4.96,Rt=0.70,Rs=1.09)', ...
          [],xy,t,mask);      % Sentinel-1
figure;imagesc(C);colorbar;


%% Test data models

% We use the current C0 as test case

C0=stmstochmodel('tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)', ...
          [],xy,t);                                          % Radarsat-2
figure;imagesc(C);colorbar;

% Pack and test

C=stmstochmodel('variance()', diag(C0) );
max(max(abs(diag(C)-diag(C0))))

C=stmstochmodel('stdev()', sqrt(diag(C0)));
max(max(abs(diag(C)-diag(C0))))

C=stmstochmodel('covmatrix(format=full)', C0);
max(max(abs(C-C0)))

stochData=C0(tril(true(size(C0))));
C=stmstochmodel('covmatrix(format=lower)', stochData);
max(max(abs(C-C0)))

[stochData,d]=symdiags(C0);
C=stmstochmodel('covmatrix(format=symdiags)', stochData);
max(max(abs(C-C0)))

%%
np=size(xy,1);
nt=size(t(:),1);
stochData=zeros(np,np,nt);
ref=cell(nt,1);
k2=0;
for k=1:nt
   k1=k2+1;
   k2=k2+np;
   stochData(1:np,1:np,k)=C0(k1:k2,k1:k2);
   ref{k}=C0(k1:k2,k1:k2);
end
Cref=blkdiag(ref{:});

C=stmstochmodel(['covmatrix(format=blkdiag,nd=' num2str(nt) ')'], stochData);
max(max(abs(C-Cref)))


%% Test recommended NAM/06-GPS model

% recommended NAM/06-GPS model continuous GNSS
C=stmstochmodel({ ...
      'iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sigma=0.0967,fs=24,ndays=3)', ...  
      'iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sigma=0.04,fs=24,ndays=3)' ...
      'spatialcov(rho=0.0887)'}, ...
      [],xy,t);   
figure;imagesc(C);colorbar

% recommended NAM/06-GPS model continuous GNSS 
C=stmstochmodel({ ...
      'iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sigma=0.0967,fs=24,ndays=1)', ...  
      'iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sigma=0.04,fs=24,ndays=1)' ...
      'spatialcov(rho=0.0887)'}, ...
      [],xy,t);   
figure;imagesc(C);colorbar

% recommended NAM/06-GPS model continuous GNSS 
C=stmstochmodel({ ...
      'iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sigma=0.0967,fs=24,ndays=1,leadup=180)', ...  
      'iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sigma=0.04,fs=24,ndays=1,leadup=180)' ...
      'spatialcov(rho=0.0887)'}, ...
      [],xy,t);   
figure;imagesc(C);colorbar

