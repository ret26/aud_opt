function hFig = plotFB(fc,betas,g_gam,fs,T,DS)

D = length(fc);
bw = audfiltbw(fc).*betas;
freqs = linspace(0,fs/2,T/2);
spec = gammatone_spectrum(fc,bw,freqs);

yspec = zeros(T,1);
yspec(floor(T/2)) = 1;
Yspec = ufilterbank(yspec,g_gam,DS);
Yspec = real(Yspec);
specY = abs(fft(Yspec));
freqs = linspace(0,fs/2,floor(T/2));

hFig = figure;
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
