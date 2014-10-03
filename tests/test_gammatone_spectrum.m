function test_suite = test_gammatone_spectrum
  initTestSuite;

% Tests: 
%
% 

%

function test_check_gradient

T = 16000;
y = zeros(T,1);
y(4000) = 1;
fs = 16000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set
DS=1; % downsampling (optional)
channels_per_erb=0.5;
filterlength = 5000; % not sure whether this is useful -- try optimising later


% Number of channels, slightly less than 1 ERB(Cambridge) per channel.
D=ceil(freqtoerb(fs/2)*channels_per_erb);

% Compute center frequencies.
%fc=erbspace(0,fs/2,D);

fc=erbspace(100,fs/2-5000,D);

% choose gammatone filter coefficients
betas = ones(1,D)*2;
g_gam=gammatonefir(fc,fs,filterlength,betas,'peakphase');

% Analysis transform
Y = ufilterbank(y,g_gam,DS);
Y = real(Y);

freqs = linspace(0,fs/2,T/2);
specY = abs(fft(Y));

bw = audfiltbw(fc).*betas;
spec = gammatone_spectrum(fc,bw,freqs);


% something odd happening at the extremes, but roughly correct
figure
hold on
for d=1:D
  plot(freqs,specY(1:T/2,d),'-k','linewidth',2)
  plot(freqs,spec(:,d),'-r')
end

%tol = 1e-4;
%assertVectorsAlmostEqual(d,0,'absolute',tol,0)


keyboard