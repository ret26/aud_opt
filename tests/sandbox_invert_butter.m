clear; 

% test inversion of butterworth filter in a stand alone file

randn('state',1);

% load sound
[y,fs] = wavread('/home/rich/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav');

% pick a short section
% y = y(8438:12190);
% y = y(8438:33630);
 y = y(8438:33630);

 
T = length(y);
Tsmooth = 100;
sinOn = sin(2*pi*[0:Tsmooth]'/(4*Tsmooth));
cosine_smooth = [sinOn;ones(T-2*Tsmooth-2,1);sinOn(end:-1:1)];
y = y.*cosine_smooth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = 5;   fc=erbspace(200,5500,D); bet = 2;
filterlength = 5000; % not sure whether this is useful -- try optimising later
betas = bet*ones(1,D); % optional widening of the filters
DS = 1;

% choose gammatone filter coefficients
g_gam=gammatonefir(fc,fs,filterlength,betas,'peakphase');

rho = 1e4;
fCutLP = 770*2/fs;
ordLP = 7;


% normalise
y = y/sqrt(var(y));

% % Analysis transform
Y = ufilterbank(y,g_gam,DS);
Y = real(Y);

[T,D] = size(Y);

% hair cell
% (soft) half wave rectify
YHW = softHWR(Y,rho);

% low pass filter
[z,p] = butter(ordLP,fCutLP);
A = zeros(T,D);
for d=1:D
  A(:,d) = filter(z,p,YHW(:,d));
end

% test to invert filter
yDelta = zeros(T,1);
yDelta(1) = 1;
aImp = filter(z,p,yDelta);

% get the FFT of the butterworth filter
spec = fft(aImp);
spec_old = abs(spec); 
max_gain = 1e3;

A2 = zeros(T,D);
for d=1:D
  A2(:,d) = ifft(fft(YHW(:,d)).*spec);
end

% show that the fft filtering method gives the same solution as
% matlab's 'filter' function

d= 5;
figure
hold on
plot(A(:,d),'-k')
plot(A2(:,d),'-r')


% test to invert filter
yDelta = zeros(2*T,1);
yDelta(1) = 1;
aImp = filter(z,p,yDelta);

% get the FFT of the butterworth filter
spec = fft(aImp);
spec_old = abs(spec); 
max_gain = 1e4;

% invert the fft based filter
spec_new = 1./abs(spec);
spec_new(spec_new>max_gain)= max_gain;


YHWrecon = zeros(T,D);
for d=1:D
  filt_chan = ifft(fft([A(:,d);A(end:-1:1,d)]).*spec_new.*exp(-i*angle(spec)));
  filt_chan = real(filt_chan(1:T));
  YHWrecon(:,d) = filt_chan;
  ind = filt_chan<0;
  YHWrecon(ind,d) = 0;
end

d = 5;
% show that the inversion is fairly good:
figure
hold on
plot(YHW(:,d),'-k')
plot(YHWrecon(:,d),'-r')

snr(YHW,YHWrecon)
