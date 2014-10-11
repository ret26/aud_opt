clear; 

% see if scale of initialisation matters and test how long we
% should optimise for
% paths that need to be added 
% pat = genpath('~/Programs/ltfat/'); addpath(pat);
% pat = genpath('~/Synchronised/'); addpath(pat);

% load sound
[y,fs] = wavread('~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav');

% pick a short section
% y = y(8438:12190);
 y = y(24480:32190);

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
hFig = plotFB(fc,betas,g_gam,fs,T,DS);
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


numIts = ones(10,1)*120;

L = 5; % number of random restarts

Ys1 = zeros(T,L);
Ys2 = zeros(T,L);

for l=1:L
  disp(['%%%%%% Progress ',num2str(l),'/',num2str(L),' %%%%%%'])
  
  yInit1 = randn(T,1)/100;
  yInit2 = yInit1*1e-5;%randn(T,1)/100;


  [ynew1,info1] = match_envelopes_V1(yInit1,ATar,g_gam,DS,fCutLP, ...
				   ordLP,rho,numIts,y,0);

  
  Ys1(:,l) = ynew1;
  objs1{l} = info1.obj;
  err_ys1{l} = info1.err_y;

  [ynew2,info2] = match_envelopes_V1(yInit2,ATar,g_gam,DS,fCutLP, ...
				     ordLP,rho,numIts,y,0);

  Ys2(:,l) = ynew2;
  objs2{l} = info2.obj;
  err_ys2{l} = info2.err_y;
  
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
  plot(Ys2(:,l),'-r')
end

subplot(2,1,2)
hold on
plot(y,'-k')

for l=1:L
  plot(Ys2(:,l),'-r')
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
