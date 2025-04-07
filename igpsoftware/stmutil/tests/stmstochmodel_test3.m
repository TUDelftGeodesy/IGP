%% Test script for stmstochmodel - iGGM comparisons
%
%   Tests using the FIGGM model with spectral index ALPHAR to right of 
%   cross-over frequency and spectral index ALPHAL to the left
%   of the cross-over frequency. The factor GAMMA (0 -> 1) defines the 
%   cross-over for AR(1).  As GAMMA tend to 1 the cross-over frequency 
%   decreases in frequency. 
%
%   A Fractionally Integrated Generalized Gauss Markov (FIGGM) model has 
%   two regions separated at some frequency, fc, which is related to the 
%   AR(1) parameter, GAMMA. To the left of fc, the spectrum has a slope in 
%   log-log space equal to ALPHAL and a slope of ALPHAR to the right of fc. 
%
%   All of the other models can be considered special cases of FIGGM. 
%
%   Model           ALPHAL     ALPHAR     GAMMA
%   -------------   -------    -------    -------
%   AR(1)/FOGM        0          -2       GAMMA
%   ARFIMA/ARFI     ALPHAL     ALPHAL-2   GAMMA
%   POWER-LAW       ALPHA      ALPHA        0       (ALPHA->ALPHAL)
%   POWER-LAW       ALPHA      ALPHA        1       (ALPHA->APLHAR)
%   GGM               0        ALPHAR     GAMMA
%   FIGGM           ALPHAL     ALPHAR     GAMMA
%
%   For instance AR(1) (also called First-Order Gauss Markov (FOGM) noise) 
%   simply has one parameter, GAMMA, which describes where the break point 
%   is; ALPHAL is 0 and ALPHAR is -2. Power-law noise which simply has one 
%   slope can arise from a number of combinations of the three parameters. 
%
%   Ref: https://www.psmsl.org/products/trends/methods.php
%
%   The autocovariance function for AR(1) decays with a decay time (also 
%   called time constant) TAU=-1/LOG(GAMMA). The cross-over frequency is
%   associated to 1/TAU = -LOG(GAMMA). [THIS IS TO BE CONFIRMED]
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

%% Test white noise models - all should give identical results

C0=stmstochmodel('WN(sigma=0.5,fs=1,ndays=1)',[],[],t);
figure;imagesc(C0);colorbar

C=stmstochmodel('FOGM(rho=Inf,sigma=0.5,fs=1,ndays=1)',[],[],t);
max(max(abs(C-C0)))

C=stmstochmodel('PL(alpha=0,sigma=0.5,fs=1,ndays=1)',[],[],t);
max(max(abs(C-C0)))

C=stmstochmodel('iGGM(alphaR=0,alphaL=0,gamma=0,sigma=0.5,fs=1,ndays=1)',[],[],t);
max(max(abs(C-C0)))

C=stmstochmodel('iGGM(alphaR=0,alphaL=0,gamma=1,sigma=0.5,fs=1,ndays=1)',[],[],t);
max(max(abs(C-C0)))


%% Test FOGM/AR(1) noise models - all should give identical results

C0=stmstochmodel('FOGM(rho=1/365.25,sigma=0.5,fs=1,ndays=1)',[],[],t);
figure;imagesc(C0);colorbar

%   AR(1) (also called First-Order Gauss Markov (FOGM) noise) simply has one 
%   parameter, GAMMA, which describes where the break point is; ALPHAL is 0 
%   and ALPHAR is -2. The autocovariance function for AR(1) decays with a 
%   decay time (also called time constant) TAU=-1/LOG(GAMMA). The cross-over 
%   frequency is associated to 1/TAU = -LOG(GAMMA), which is equal to the 
%   parameter RHO in the FOGM model. The SIGMA must also be adjusted,
%
%      GAMMA = EXP(-1/TAU) 
%      SIGMA = SIGMA0 / SQRT(TAU/2)
% 
%  with RHO=1/TAU and SIGMA0 the parameters for the FOGM.

gamma=exp(-1/365.25);
num2str(gamma,8)

%C=stmstochmodel('iGGM(alphaR=-2,alphaL=0,gamma=0.997266,sigma=0.5/sqrt(365.25/2),fs=1,ndays=1,leadup=1000)',[],[],t);
C=stmstochmodel(['iGGM(alphaR=-2,alphaL=0,gamma=' num2str(gamma,8) ',sigma=0.5/sqrt(365.25/2),fs=1,ndays=1,leadup=1000)'],[],[],t);
figure;imagesc(C);colorbar
max(max(abs(C-C0)))


%% Test power law models - all should give identical results

alpha=-.8;
C0=stmstochmodel(['PL(alpha=' num2str(alpha) ',sigma=0.5,fs=1,ndays=1)'],[],[],t);
figure;imagesc(C0);colorbar

C=stmstochmodel(['iGGM(alphaR=' num2str(alpha) ',alphaL=' num2str(alpha) ',gamma=1,sigma=0.5,fs=1,ndays=1)'],[],[],t);
max(max(abs(C-C0)))

C=stmstochmodel(['iGGM(alphaR=' num2str(alpha) ',alphaL=' num2str(alpha) ',gamma=0,sigma=0.5,fs=1,ndays=1)'],[],[],t);
max(max(abs(C-C0)))

%%  CATS scaling 
%
%     sigma = dt^(-alpha/4) * sigCATS
%
alpha=-.8281
sigma=0.76
dt=1/365.25;

dt^(-alpha/4)

sigCATS=sigma*dt^(alpha/4)
