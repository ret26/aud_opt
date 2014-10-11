function test_suite = test_invert_butter
  initTestSuite;

% Tests: 
%
% 
% [obj,dobj] = getObjV1(y,Atar,g_gam,DS,fCutLP,ordLP,rho)
%

function test_inversion

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


max_gain = 1e4;
YHWrecon = invert_LPF(A,rho,fCutLP,ordLP,max_gain);


d = 1;
% show that the inversion is fairly good:
figure
subplot(1,3,1)
hold on
plot(YHW(:,d),'-k')
plot(YHWrecon(:,d),'-r')

d = 3;
subplot(1,3,2)
hold on
plot(YHW(:,d),'-k')
plot(YHWrecon(:,d),'-r')

d = 5;
subplot(1,3,3)
hold on
plot(YHW(:,d),'-k')
plot(YHWrecon(:,d),'-r')

snrs = snr(YHW,YHWrecon);


%tol = 1e-4;
%assertVectorsAlmostEqual(d,0,'absolute',tol,0)


