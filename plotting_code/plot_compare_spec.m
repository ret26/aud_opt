function figh = plot_compare_spec(ynew,y,fs,g_gam,DS)

T = length(y);
D = length(g_gam);

figh = figure;
hold on
disc = y - ynew;
spec_disc = abs(fft(disc));
spec_y = abs(fft(y));
spec_ynew = abs(fft(ynew));


yspec = zeros(T,1);
yspec(floor(T/2)) = 1;
Yspec = ufilterbank(yspec,g_gam,DS);
Yspec = real(Yspec);
specY = abs(fft(Yspec));
freqs = linspace(0,fs/2,floor(T/2));
for d=1:D
  plot(freqs,specY(1:floor(T/2),d),'-','linewidth',2,'color',[1,1,1]*0.8)
end

scale = max([spec_y;spec_ynew]);

freq = linspace(0,fs/2,floor(T/2));
plot(freq,spec_y(1:floor(T/2))/scale,'-k')
plot(freq,spec_ynew(1:floor(T/2))/scale,'-r')
plot(freq,spec_disc(1:floor(T/2))/scale)

ylabel('spectrum')
xlabel('frequency')
set(gca,'yscale','log')