clear; 

% compares reconstruction from the envelopes to reconstruction from
% the reconstructed half wave rectified filter coefficients

randn('state',1);

% load sound
[y,fs] = wavread('~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav');

% pick a short section
% y = y(8438:12190);
% y = y(8438+3887:8640+8438);

% pick a long section
 y = y(8438:33630);

% normalise
y = y/sqrt(var(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set auditory model 
DS = 1; % downsampling (optional)

T=length(y);
% Number of channels, slightly less than 1 ERB(Cambridge) per channel.

%channels_per_erb = 0.2; D=ceil(freqtoerb(fs/2)*channels_per_erb); %fc=erbspace(50,fs/2-1000,D);

%D = 20;   fc=erbspace(200,5500,D); bet = 1;
%D = 10;   fc=erbspace(100,5500,D); bet = 1;
D = 5;   fc=erbspace(100,5500,D); bet = 2;
%D = 3;   fc=erbspace(100,5500,D); bet = 3;
%D = 4;   fc=erbspace(100,5500,D); bet = 2;


filterlength = 5000; % not sure whether this is useful -- try optimising later
betas = bet*ones(1,D); % optional widening of the filters

% choose gammatone filter coefficients
g_gam=gammatonefir(fc,fs,filterlength,betas,'peakphase');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot filter bank
hFig = plotFB(fc,betas,g_gam,fs,T,DS);
drawnow;
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hair cell
% (soft) half wave rectify
rho = 1e4;
fCutLP = 770*2/fs;
ordLP = 7;

% compute targets 
[ATar,YHWTar,YTar] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho);

max_gain = 1e3;
YHWTar = invert_LPF(ATar,rho,fCutLP,ordLP,max_gain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimise new signal

ynew = 0*randn(T,1);
numIts = [1,1,1,1,1,1,2,2,2,2,16,32,64,128,128];

numFrames = length(numIts);
Ys = zeros(T,numFrames+2);
Ys(:,1) = ynew;

for fr = 1:numFrames
  
  numItsCur = numIts(fr);
  
  [ynew,infoInit2] = match_HWR_filt_V1(ynew,YHWTar,g_gam,DS,rho, ...
					 numItsCur,y);
  Ys(:,fr+1) = ynew;

  figure
  hold on;
  plot(y,'-k')
  plot(ynew,'-r')
  drawnow;
  pause(0.1)
  
end

numIts = 160;
[ynew,info2] = match_envelopes_V1(ynew,ATar,g_gam,DS,fCutLP,ordLP,rho,numIts,y);

Ys(:,numFrames+2) = ynew;

figure
hold on;
plot(y,'-k')
plot(ynew,'-r')

drawnow;
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolate to add more frames

[T,K] = size(Ys);

Ysnew = zeros(T,2*K-1);

for k=1:K
    kk=2*(k-1)+1;
    Ysnew(:,kk) = Ys(:,k);
    if k~=K
      Ysnew(:,kk+1) = (Ys(:,k) + Ys(:,k+1))/2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveDir = '~/Data/aud_opt/mk_video_V1/vid3/';
printOutput=1;
saveName = 'first';

DSy = 2; % amount to downsample y
DSA = 100; % amount to downsample A

for fr=1:size(Ysnew,2)
  plot_video_frame_V1(y,Ysnew,g_gam,DS,fCutLP,ordLP,rho,fc,fs,fr, ...
		      DSy,DSA,saveDir,saveName,printOutput)
  close all;
end


printOutput=1;
saveName = 'close_up';
%ind = [8438+3887:8640+8438];
ind = [1+200:12190-8438-200];

DSy = 1; % amount to downsample y
DSA = 5; % amount to downsample A


for fr=1:size(Ysnew,2)
  plot_video_frame_V1(y(ind),Ysnew(ind,:),g_gam,DS,fCutLP,ordLP, ...
		      rho,fc,fs,fr,DSy,DSA,saveDir,saveName,printOutput)
  
  close all;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save sounds
savefile = '~/Data/aud_opt/mk_video_V1/';
savebasename = ['video_',num2str(bet),'_D_',num2str(D),'_']
save([savefile,savebasename,'.mat']);
