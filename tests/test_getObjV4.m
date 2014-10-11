function test_suite = test_getObjV4
  initTestSuite;

% Tests: 
%
% 
% function [obj,dobj] = getObjV4(y,Atar,g_gam,DS,fCutLP,ordLP,rho,kappa1,kappa2,yorig)
%

function test_check_gradient

T = 300;
y = randn(T,1);
yorig = randn(T,1);
fs = 1000;

% filterbank parameters 
DS=1; % downsampling (optional)
channels_per_erb=0.25; % normally 2
filterlength = 208; % not sure whether this is useful -- try optimising later

% hair cell paramters
rho = 1e4;
fCutLP = 77*2/fs;
ordLP = 7;

% prior strength
kappa1 = rand;
kappa2 = rand;

% Number of channels, slightly less than 1 ERB(Cambridge) per channel.
D=ceil(freqtoerb(fs/2)*channels_per_erb);
fc=erbspace(0,fs/2,D); % Compute center frequencies.
g_gam=gammatonefir(fc,fs,filterlength,'peakphase'); % choose gammatone filter coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Atar = rand(T,D);

[obj,dobj] = getObjV4(y,Atar,g_gam,DS,fCutLP,ordLP,rho,kappa1,kappa2,yorig);

delta = 1e-5;
yinit = randn(T,1);

d=checkgrad('getObjV4',yinit,delta,Atar,g_gam,DS,fCutLP,ordLP,rho,kappa1,kappa2,yorig);
d
tol = 1e-4;
assertVectorsAlmostEqual(d,0,'absolute',tol,0)


