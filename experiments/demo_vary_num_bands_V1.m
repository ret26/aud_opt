clear; 

% alter the number of bands and measure the error w.r.t. y

% load sound
[y,fs] = wavread('/home/rich/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav');

% pick a short section
 y = y(8438:12190);

% pick a long section
% y = y(8438:33630);

% normalise
y = y /sqrt(var(y));

T=length(y);

DS=1; % downsampling (optional)

Ds = [1,3,6,9];

% hair cell
% (soft) half wave rectify
rho = 1e4;
fCutLP = 770*2/fs;
ordLP = 7;
filterlength = 5000; % not sure whether this is useful -- try optimising later

% numbers of iterations and restarts
numIts = ones(15,1)*40;
L = 5; % number of random restarts
Ys = zeros(T,L,M);
err_ys = zeros(length(numIts,L,M);


for m=1:M
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set auditory model 

  % Number of channels, slightly less than 1 ERB(Cambridge) per channel.
  D = Ds(m); %=ceil(freqtoerb(fmax)*channels_per_erb(m));
  
  % Compute center frequencies.
  fc=erbspace(50,fmax,D);
  betas = bet(m)*ones(1,D); % optional widening of the filters 
  
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

  [ATar,YHWTar,YTar] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for l=1:L
    disp(['%%%%%% Progress ',num2str(l),'/',num2str(L),' %%%%%%'])
    
    yInit = randn(T,1);

    [ynew,info] = match_envelopes_V1(yInit,ATar,g_gam,DS,fCutLP, ...
				       ordLP,rho,numIts,y,0);
   
    
    Ys(:,l,m) = ynew;
    objs{l,m} = info.obj;
    err_ys(:,l,m) = info.err_y;
    
  end
end

% [A,YHW,Y] = aud_mod_V1(ynew,g_gam,DS,fCutLP,ordLP,rho);
figure
subplot(1,2,1)
hold on
ylabel('objective')

for l=1:L
  plot(objs1{l},'-k')
  plot(objs2{l},'-r')
end

set(gca,'yscale','log')
legend('condition 1','condition 2')

subplot(1,2,2)
hold on

for l=1:L
  plot(err_ys1{l},'-k')
  plot(err_ys2{l},'-r')
end
ylabel('error y')
set(gca,'yscale','log')




%K = 20;

% figure
% subplot(2,1,1)

% tol = 1e-5;
% ATarPlot = ATar;
% ATarPlot(ATar<tol) = tol;
% surf(log(ATarPlot(1:K:end,:))','edgecolor','none')
% view(0,90)
% set(gca,'xlim',[1,T/K],'ylim',[1,D])

% subplot(2,1,2)
% APlot = A;
% APlot(A<tol) = tol;
% surf(log(APlot(1:K:end,:))','edgecolor','none')
% view(0,90)
% set(gca,'xlim',[1,T/K],'ylim',[1,D])


figure
subplot(2,1,1)
hold on
plot(y,'-k')

for l=1:L
  plot(Ys1{l},'-r')
end

subplot(2,1,2)
hold on
plot(y,'-k')

for l=1:L
  plot(Ys1{l},'-r')
end

mn_err_y1 = 0;
mn_err_y2 = 0;
mn_objs1 = 0;
mn_objs2 = 0;

for l=1:L
  mn_err_y1 = mn_err_y1+1/L*err_ys1{1}(end);
  mn_err_y2 = mn_err_y2+1/L*err_ys2{1}(end);
  mn_objs1 = mn_objs1+1/L*objs1{1}(end);
  mn_objs2 = mn_objs2+1/L*objs2{1}(end); 
end

disp('mean objectives reached')
[mn_objs1,mn_objs2]

disp('mean squared discrepancy from true signal')
[mn_err_y1, mn_err_y2]

disp('variance in resulting solutions')
[mean(var(Ys1')),mean(var(Ys2'))]

%yscale = max(abs([y;ynew;yvoc]));
%savefile = '/home/rich/Synchronised/aud_opt/sounds/';
%wavwrite(0.95*y/yscale,fs,[savefile,'original.wav'])
%wavwrite(0.95*ynew/yscale,fs,[savefile,'chimera1.wav'])
%wavwrite(0.95*yvoc/yscale,fs,[savefile,'chimera2.wav'])
