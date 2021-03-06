clear; 

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
%D = 10;   fc=erbspace(200,5500,D); bet = 1;
%D = 3;   fc=erbspace(200,5500,D); bet = 3;
%D = 4;   fc=erbspace(200,5500,D); bet = 2;
%D = 5;   fc=erbspace(100,5500,D); bet = 2;
D = 3;   fc=erbspace(100,5500,D); bet = 3;

%D = 3; fc = [100,500,4000]; bet = 3;
%D = 2; fc = [200,800]; bet = 5;
%D = 1; fc = [100]; bet = 7;

filterlength = 5000; % not sure whether this is useful -- try optimising later
betas = bet*ones(1,D); % optional widening of the filters

% choose gammatone filter coefficients
g_gam=gammatonefir(fc,fs,filterlength,betas,'peakphase');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot filter bank
hFig = plotFB(fc,betas,g_gam,fs,T,DS);

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

max_gain = 1e4;
YHWTar = invert_LPF(ATar,rho,fCutLP,ordLP,max_gain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimise new signal

yInit = randn(T,1)*0;
numIts = ones(10,1)*20;
%numIts = ones(40,1)*40;

[ynew,info] = match_HWR_filt_V1(yInit,YHWTar,g_gam,DS,rho,numIts,y,fCutLP,ordLP,ATar);
%[ynew,info] = match_envelopes_V1(ynew,ATar,g_gam,DS,fCutLP,ordLP,rho,numIts,y);

[A,YHW,Y] = aud_mod_V1(ynew,g_gam,DS,fCutLP,ordLP,rho);

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

[snr_y,snr_A] = compute_sig_sim(y,ynew,fs);
[snr_yvoc1,snr_Avoc1] = compute_sig_sim(y,yvoc1,fs);
[snr_yvoc2,snr_Avoc2] = compute_sig_sim(y,yvoc2,fs);
[snr_yvoc3,snr_Avoc3] = compute_sig_sim(y,yvoc3,fs);

snr_Atrain = snr(ATar,A);

disp(['SNR of optimised reconstruction ',num2str(snr_y),'dB'])
disp(['SNR of vocoded reconstruction 1 ',num2str(snr_yvoc1),'dB'])
disp(['SNR of vocoded reconstruction 2 ',num2str(snr_yvoc2),'dB'])
disp(['SNR of vocoded reconstruction 3 ',num2str(snr_yvoc3),'dB'])

disp(['SNR of optimised reconstruction envelopes ',num2str(mean(snr_A)),'dB'])
disp(['SNR of vocoded recon 1 envelopes ',num2str(mean(snr_Avoc1)),'dB'])
disp(['SNR of vocoded recon 2 envelopes ',num2str(mean(snr_Avoc2)),'dB'])
disp(['SNR of vocoded recon 2 envelopes ',num2str(mean(snr_Avoc3)),'dB'])
disp(['SNR of optimised reconstruction on training enveopes ',num2str(mean(snr_Atrain)),'dB'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure
subplot(1,2,1)
hold on
ylabel('objective')
plot(info.obj)
set(gca,'yscale','log')
objvoc1 = getObjV1(yvoc1,ATar,g_gam,DS,fCutLP,ordLP,rho);
objvoc2 = getObjV1(yvoc2,ATar,g_gam,DS,fCutLP,ordLP,rho);
objvoc3 = getObjV1(yvoc3,ATar,g_gam,DS,fCutLP,ordLP,rho);
objnoi = getObjV1(randn(T,1),ATar,g_gam,DS,fCutLP,ordLP,rho);
plot(length(info.obj),objvoc1,'ok')
plot(length(info.obj),objvoc2,'ob')
plot(length(info.obj),objvoc3,'om')
plot(length(info.obj),objnoi,'or')
legend('optimised','noise voc. 1','noise voc. 2','noise voc. 3','noise','location','northwest')

subplot(1,2,2)
hold on
plot(info.err_y)
ylabel('error y')
set(gca,'yscale','log')
objyvoc1 = mean((y-yvoc1).^2);
objyvoc2 = mean((y-yvoc2).^2);
objyvoc3 = mean((y-yvoc3).^2);
objynoi = mean((y-randn(T,1)).^2);
%plot(length(info.err_y),objyvoc,'ok')
%plot(length(info.err_y),objynoi,'or')

figh = plot_snr_channels(ynew,y,fc,fs,g_gam,DS,fCutLP,ordLP,rho);

figh = plot_snr_channels_voc(ynew,yvoc2,y,fs);

figh = plot_compare_spec(ynew,y,fs,g_gam,DS);

figh = plot_disc(ynew,y,fs);

figh = plot_compare_sigs(ynew,y,fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save sounds
yscale = max(abs([y;ynew;yvoc1;yvoc2;yvoc3]));
savefile = '~/data/aud_opt/demo_short_speech_V2/';
savebasename = ['sentence_bet_',num2str(bet),'_D_',num2str(D),'_']
wavwrite(0.95*y/yscale,fs,[savefile,savebasename,'original.wav'])
wavwrite(0.95*ynew/yscale,fs,[savefile,savebasename,'chimera1.wav'])
wavwrite(0.95*yvoc1/yscale,fs,[savefile,savebasename,'chimera2.wav'])
wavwrite(0.95*yvoc2/yscale,fs,[savefile,savebasename,'chimera3.wav'])
wavwrite(0.95*yvoc3/yscale,fs,[savefile,savebasename,'chimera4.wav'])
