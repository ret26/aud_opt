function [snr_y,snr_A] = compute_sig_sim(y,ynew,fs)

% function [snr_y,snr_A] = compute_sig_sim(y,ynew)

% SNR waveform
snr_y = snr(y,ynew);

% SNR full filerbank

channels_per_erb = 1;
filterlength = 5000; % not sure whether this is useful -- try optimising later

T=length(y);
% Number of channels, slightly less than 1 ERB(Cambridge) per channel.
D=ceil(freqtoerb(fs/2)*channels_per_erb);

% Compute center frequencies.
fc=erbspace(50,fs/2-1000,D);

betas = 1*ones(1,D); % optional widening of the filters

% choose gammatone filter coefficients
g_gam=gammatonefir(fc,fs,filterlength,betas,'peakphase');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot filter bank
% hFig = plotFB(fc,betas,g_gam,fs,T,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keyboard

% hair cell
% (soft) half wave rectify
rho = 1e4;
fCutLP = 770*2/fs;
ordLP = 7;

[ATar,YHWTar,YTar] = aud_mod_V1(y,g_gam,1,fCutLP,ordLP,rho);
[A,YHW,Y] = aud_mod_V1(ynew,g_gam,1,fCutLP,ordLP,rho);

snr_A = snr(ATar,A);
