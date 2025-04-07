%% Script to test robust moving average filter rmafilt.m

%% Generate test dataset with steps and outlies

% Set Length of the series, simulated steps and outliers 

n=1000;
sigma0=1.5;
            %  epoch magnitude
simsteps    = [ 300     3 ; ...
                550    -5 ; ...
                610     3 ; ...
                630    -5 ];
simoutliers = [ 200     5 ; ...
                500     6 ; ...
                600    -7 ];

% Simulate the test dataset

y=sigma0*randn(n,1);
y(simoutliers(:,1))=simoutliers(:,2);
for k=1:size(simsteps,1)
  s=zeros(size(y));
  s(simsteps(k,1):end,1)=1;
  y=y+simsteps(k,2)*s;
end

%% Call rmafilt

[ym,sy,idxoutlier,idxstep,steps]=rmafilt(y, 11, 'verbose', 2, 'simoutliers', simoutliers ,'simsteps', simsteps );
