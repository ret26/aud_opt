clear;

% demo of the simple version of the auditory model
%pat = genpath('~/Programs/ltfat/'); addpath(pat);
%pat = genpath('~/Synchronised/'); addpath(pat);

% load sound
[y,fs] = wavread('/home/rich/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav');

y = y(7000:2.1*fs);

DS=1; % downsampling (optional)
channels_per_erb=2;
filterlength = 5000; % not sure whether this is useful -- try optimising later

% Determine minimal transform length
T=length(y);
L=ceil(filterlength/DS)*DS;

% Number of channels, slightly less than 1 ERB(Cambridge) per channel.
D=ceil(freqtoerb(fs/2)*channels_per_erb);

% Compute center frequencies.
fc=erbspace(0,fs/2,D);

% choose gammatone filter coefficients
g_gam=gammatonefir(fc,fs,filterlength,'peakphase');

% apply the simple auditory model
[A,Yhw,Y] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho);

% check the filtering works
figure;
d=D;
hold on
plot(real(Y(:,d)),'-k')
plot(Yhw(:,d),'-b')
plot(A(:,d),'-r')

spec1 = abs(fft(A(:,d)));
spec2 = abs(fft(Yhw(:,d)));

figure
hold on
freqs = linspace(0,fs/2,floor(T/2));

plot(freqs,log(spec2(1:floor(T/2))+1/T))
plot(freqs,log(spec1(1:floor(T/2))+1/T),'-r')


% dynrange_for_plotting=50; 
% figure(1);
% plotfilterbank(Y,DS,fc,fs,dynrange_for_plotting,'audtick');
