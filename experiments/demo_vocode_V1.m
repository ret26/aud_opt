clear; 

% compares to three different vocoding methods:

% chimera 1: filter white noise, impose envelopes, sum to produce signal
%
% chimera 2: for each band, combine envelope and white noise, pass
%            through band, then sum
%
% chimera 3: do both of the above i.e. combine envelopes with filtered
% noise and then, for each band, put the result through the
% corresponding filter
%
% method 2 seems to be the best one (and this is the one used in
% the previous papers)

randn('state',1);

% load sound
sound1 = '/home/rich/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav';
[y,fs] = wavread(sound1);

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

%channels_per_erb = 0.2;
D = 3; fc = [100,500,4000];
% Compute center frequencies.
D= 20; fc=erbspace(50,fs/2-1000,D);


filterlength = 5000; % not sure whether this is useful -- try optimising later

bet = 2;
betas = bet*ones(1,D); % optional widening of the filters

% choose gammatone filter coefficients
g_gam=gammatonefir(fc,fs,filterlength,betas,'peakphase');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot filter bank
hFig = plotFB(fc,betas,g_gam,fs,T,DS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keyboard

% hair cell
% (soft) half wave rectify
rho = 1e4;
fCutLP = 770*2/fs;
ordLP = 7;

% compute targets 
[ATar,YHWTar,YTar] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make noise vocoded chimera

% chimera 1: filter white noise, impose envelopes, sum to produce signal
[ANoi,YHWNoi,YNoi] = aud_mod_V1(randn(T,1),g_gam,DS,fCutLP,ordLP,rho);

Yvoc1 = YNoi.*ATar;
yvoc1 = sum(Yvoc1,2);
yvoc1 = yvoc1/sqrt(var(yvoc1));
[Avoc1,YHWvoc1,Yvoc1] = aud_mod_V1(yvoc1,g_gam,DS,fCutLP,ordLP,rho);

% chimera 2: for each band, combine envelope and white noise, pass
%            through band, then sum

yvoc2 = zeros(T,1);
for d=1:D
  yCur = randn(T,1).*ATar(:,d);
  yCur = ufilterbank(yCur,g_gam(d),DS);
  yCur = real(yCur);
  yvoc2 = yvoc2+yCur; 
end

yvoc2 = yvoc2/sqrt(var(yvoc2));
[Avoc2,YHWvoc2,Yvoc2] = aud_mod_V1(yvoc2,g_gam,DS,fCutLP,ordLP,rho);

% chimera 3: do both of the above i.e. combine envelopes with filtered
% noise and then, for each band, put the result through the
% corresponding filter

[ANoi,YHWNoi,YNoi] = aud_mod_V1(randn(T,1),g_gam,DS,fCutLP,ordLP,rho);

Yvoc3 = YNoi.*ATar;

yvoc3 = zeros(T,1);
for d=1:D
  yCur = ufilterbank(Yvoc3(:,d),g_gam(d),DS);
  yCur = real(yCur);
  yvoc3 = yvoc3+yCur; 
end

yvoc3 = yvoc3/sqrt(var(yvoc3));
[Avoc3,YHWvoc3,Yvoc3] = aud_mod_V1(yvoc3,g_gam,DS,fCutLP,ordLP,rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute similarities to original signal

[snr_yvoc1,snr_Avoc1] = compute_sig_sim(y,yvoc1,fs);
[snr_yvoc2,snr_Avoc2] = compute_sig_sim(y,yvoc2,fs);
[snr_yvoc3,snr_Avoc3] = compute_sig_sim(y,yvoc3,fs);

snr_Avoc1_train = snr(ATar,Avoc1);
snr_Avoc2_train = snr(ATar,Avoc2);
snr_Avoc3_train = snr(ATar,Avoc3);

disp(['SNR of vocoded recon 1 train envelopes ',num2str(mean(snr_Avoc1_train)),'dB'])
disp(['SNR of vocoded recon 2 train envelopes ',num2str(mean(snr_Avoc2_train)),'dB'])
disp(['SNR of vocoded recon 3 train envelopes ',num2str(mean(snr_Avoc3_train)),'dB'])

disp(['SNR of vocoded reconstruction 1 ',num2str(snr_yvoc1),'dB'])
disp(['SNR of vocoded reconstruction 2 ',num2str(snr_yvoc2),'dB'])
disp(['SNR of vocoded reconstruction 3 ',num2str(snr_yvoc3),'dB'])

disp(['SNR of vocoded recon 1 all envelopes ',num2str(mean(snr_Avoc1)),'dB'])
disp(['SNR of vocoded recon 2 all envelopes ',num2str(mean(snr_Avoc2)),'dB'])
disp(['SNR of vocoded recon 3 all envelopes ',num2str(mean(snr_Avoc3)),'dB'])


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % save sounds
% yscale = max(abs([y;ynew;yvoc1;yvoc2;yvoc3]));
% savefile = '/home/rich/Data/aud_opt/demo_short_speech_V1/';
% savebasename = ['sentence_bet_',num2str(bet),'_D_',num2str(D),'_']
% wavwrite(0.95*y/yscale,fs,[savefile,savebasename,'original.wav'])
% wavwrite(0.95*ynew/yscale,fs,[savefile,savebasename,'chimera1.wav'])
% wavwrite(0.95*yvoc1/yscale,fs,[savefile,savebasename,'chimera2.wav'])
% wavwrite(0.95*yvoc2/yscale,fs,[savefile,savebasename,'chimera3.wav'])
% wavwrite(0.95*yvoc3/yscale,fs,[savefile,savebasename,'chimera4.wav'])
