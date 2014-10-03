clear; 

% load sound
[y,fs] = wavread('/home/rich/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav');

% pick a short section
 y = y(8438:12190);

% pick a long section
% y = y(8438:33630);

% normalise
y = y /sqrt(var(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set auditory model 
DS=1; % downsampling (optional)
channels_per_erb=0.2;
filterlength = 5000; % not sure whether this is useful -- try optimising later

T=length(y);
% Number of channels, slightly less than 1 ERB(Cambridge) per channel.
D=ceil(freqtoerb(fs/2)*channels_per_erb);
% Compute center frequencies.
fc=erbspace(50,fs/2-1000,D);

betas = 4*ones(1,D); % optional widening of the filters 

% choose gammatone filter coefficients
g_gam=gammatonefir(fc,fs,filterlength,betas,'peakphase');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot filter bank
bw = audfiltbw(fc).*betas;
freqs = linspace(0,fs/2,T/2);
spec = gammatone_spectrum(fc,bw,freqs);

yspec = zeros(T,1);
yspec(floor(T/2)) = 1;
Yspec = ufilterbank(yspec,g_gam,DS);
Yspec = real(Yspec);
specY = abs(fft(Yspec));
freqs = linspace(0,fs/2,floor(T/2));

figure
subplot(2,1,1)
hold on
for d=1:D
  plot(freqs,specY(1:floor(T/2),d),'-r','linewidth',2)
  plot(freqs,spec(:,d),'-k')
end

subplot(2,1,2)
hold on
for d=1:D
  plot(freqs,specY(1:floor(T/2),d),'-k','linewidth',2)
  plot(freqs,spec(:,d),'-k')
end

set(gca,'xscale','log','xlim',[10,fs/2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keyboard

% hair cell
% (soft) half wave rectify
rho = 1e4;
fCutLP = 770*2/fs;
ordLP = 7;

[ATar,YHWTar,YTar] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho);
%imagesc(log(ATar-min(ATar(:))+1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yInit = randn(T,1);
numIts = ones(6,1)*40;

[ynew,info] = match_envelopes_V1(yInit,ATar,g_gam,DS,fCutLP,ordLP,rho,numIts,y);

[A,YHW,Y] = aud_mod_V1(ynew,g_gam,DS,fCutLP,ordLP,rho);

% make noise vocoded chimera
[ANoi,YHWNoi,YNoi] = aud_mod_V1(randn(T,1),g_gam,DS,fCutLP,ordLP,rho);
yvoc = sum(YNoi.*ATar,2);
yvoc = yvoc/sqrt(var(yvoc));

figure
subplot(1,2,1)
hold on
ylabel('objective')
plot(info.obj)
set(gca,'yscale','log')
objvoc = getObjV1(yvoc,ATar,g_gam,DS,fCutLP,ordLP,rho);
objnoi = getObjV1(randn(T,1),ATar,g_gam,DS,fCutLP,ordLP,rho);
plot(length(info.obj),objvoc,'ok')
plot(length(info.obj),objnoi,'or')
legend('optimised','noise vocoded','noise')

subplot(1,2,2)
hold on
plot(info.err_y)
ylabel('error y')
set(gca,'yscale','log')
objyvoc = mean((y-yvoc).^2);
objynoi = mean((y-randn(T,1)).^2);
plot(length(info.err_y),objyvoc,'ok')
plot(length(info.err_y),objynoi,'or')

K = 20;
figure
subplot(2,1,1)

tol = 1e-5;
ATarPlot = ATar;
ATarPlot(ATar<tol) = tol;
surf(log(ATarPlot(1:K:end,:))','edgecolor','none')
view(0,90)
set(gca,'xlim',[1,T/K],'ylim',[1,D])

subplot(2,1,2)
APlot = A;
APlot(A<tol) = tol;
surf(log(APlot(1:K:end,:))','edgecolor','none')
view(0,90)
set(gca,'xlim',[1,T/K],'ylim',[1,D])


figure
hold on
plot(y,'-k')
plot(ynew,'-r')

yscale = max(abs([y;ynew;yvoc]));
savefile = '/home/rich/Synchronised/aud_opt/sounds/';
wavwrite(0.95*y/yscale,fs,[savefile,'original.wav'])
wavwrite(0.95*ynew/yscale,fs,[savefile,'chimera1.wav'])
wavwrite(0.95*yvoc/yscale,fs,[savefile,'chimera2.wav'])
