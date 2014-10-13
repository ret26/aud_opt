clear; 

randn('state',1);

% load sound
[y,fs] = wavread('~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav');

% pick a short section
% y = y(8438:12190);
 y = y(8438+3887:8640+8438);

% pick a long section
% y = y(8438:33630);

% normalise
y = y/sqrt(var(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set auditory model 
DS = 1; % downsampling (optional)

T=length(y);
% Number of channels, slightly less than 1 ERB(Cambridge) per channel.

%channels_per_erb = 0.2; D=ceil(freqtoerb(fs/2)*channels_per_erb); %fc=erbspace(50,fs/2-1000,D);

%D = 20;   fc=erbspace(200,5500,D); bet = 1;
%D = 10;   fc=erbspace(200,5500,D); bet = 1;
%D = 3;   fc=erbspace(200,5500,D); bet = 3;
%D = 4;   fc=erbspace(200,5500,D); bet = 2;
D = 5;   fc=erbspace(200,2500,D); bet = 2;
%D = 3;   fc=erbspace(100,5500,D); bet = 3;

%D = 3; fc = [100,500,4000]; bet = 3;
%D = 2; fc = [200,800]; bet = 5;
%D = 1; fc = [100]; bet = 7;

filterlength = 5000; % not sure whether this is useful -- try optimising later
betas = bet*ones(1,D); % optional widening of the filters

% choose gammatone filter coefficients
g_gam=gammatonefir(fc,fs,filterlength,betas,'peakphase');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow;
pause(0.1);

% hair cell
% (soft) half wave rectify
rho = 1e4;
fCutLP = 770*2/fs;
ordLP = 7;

% compute targets 
[ATar,YHWTar,YTar] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho);


scale = [1,1,4,6,3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot filter bank
hFig1 = plotFB(fc,betas,g_gam,fs,T,DS);


hFig2 =figure;
plot(y,'-k')
set(gca,'visible','off')

sep = 2*max(abs(YTar(:)));

hFig3 =figure;
hold on
for d=1:D
  plot(YTar(:,d)*scale(d)+sep*(d-1),'-k','linewidth',2);
end
set(gca,'visible','off')

hFig4 =figure;
hold on
for d=1:D
  plot(YHWTar(:,d)*scale(d)+sep*(d-1),'-k','linewidth',1);
end
set(gca,'visible','off')


hFig5 =figure;
hold on
for d=1:D
  plot(ATar(:,d)*scale(d)+sep*(d-1),'-k','linewidth',1);
end
set(gca,'visible','off')

