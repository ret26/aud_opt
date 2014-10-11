clear; 
randn('state',1);

% alter the number of bands and measure the error w.r.t. y

% load sound
soundpath1 = '/home/rich/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav';

soundpath2 = '~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav';

savedir = '~/data/aud_opt/demo_vary_num_bands_V1/'
savename = 'first_effort';

[y,fs] = wavread(soundpath2);

% pick a short section
% y = y(8438:12190);

% pick a long section
 y = y(8438:33630);

% normalise
y = y /sqrt(var(y));

T=length(y);

DS=1; % downsampling (optional)

Ds = [2,3,4,5,6,10,20];
bet = [3,3,2.5,2.5,2,1,1];

%fmin = 100;
%fmax = 6000;

fmin = 200;
fmax = 5500;

% hair cell
% (soft) half wave rectify
rho = 1e4;
fCutLP = 770*2/fs;
ordLP = 7;
filterlength = 5000; % not sure whether this is useful -- try optimising later

% numbers of iterations and restarts
numIts = ones(20,1)*80;
L = 1; % number of random restarts

M = length(Ds); % number of different numbers of channels to use

Ys = zeros(T,L,M);
Ysvoc = zeros(T,L,M);
err_ys = zeros(length(numIts),L,M);

snr_y = zeros(L,M);
snr_A = zeros(L,M);
snr_yvoc = zeros(L,M);
snr_Avoc  = zeros(L,M);
snr_Atrain = zeros(L,M);

for m=1:M
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set auditory model 

  % Number of channels, slightly less than 1 ERB(Cambridge) per channel.
  D = Ds(m); %=ceil(freqtoerb(fmax)*channels_per_erb(m));
  
  % Compute center frequencies.
  fc=erbspace(fmin,fmax,D);
  betas = bet(m)*ones(1,D); % optional widening of the filters 
  
  % choose gammatone filter coefficients
  g_gam=gammatonefir(fc,fs,filterlength,betas,'peakphase');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot filter bank
  %hFig = plotFB(fc,betas,g_gam,fs,T,DS);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %keyboard

  [ATar,YHWTar,YTar] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for l=1:L
    disp(['%%%%%% Progress ',num2str(l),'/',num2str(L),' %%%%%%'])
    
    yInit = 0*randn(T,1);

    [ynew,info] = match_envelopes_V1(yInit,ATar,g_gam,DS,fCutLP, ...
				       ordLP,rho,numIts,y,0);
   
    
    Ys(:,l,m) = ynew;
    objs{l,m} = info.obj;
    err_ys(:,l,m) = info.err_y;
    
    yvoc = noise_vocode(ATar,g_gam,DS,fCutLP,ordLP,rho);
    Ysvoc(:,l,m) = yvoc;

    info.obj(end)
    
    % various metrics
    
    [A,YHW,Y] = aud_mod_V1(ynew,g_gam,DS,fCutLP,ordLP,rho);
    snr_Atraincur = snr(ATar,A);
    snr_Atrain(l,m) = mean(snr_Atraincur);
    
    [snr_y(l,m),snr_Acur] = compute_sig_sim(y,ynew,fs);
    snr_A(l,m) = mean(snr_Acur);
    
    [snr_yvoc(l,m),snr_Avoccur] = compute_sig_sim(y,yvoc,fs);
    snr_Avoc(l,m) =  mean(snr_Avoccur);
    
    disp('snr y')
    round(snr_y*10)/10
    
    disp('snr A (all)')
    round(snr_A*10)/10
    
    disp('snr A (target)')
    round(snr_Atrain*10)/10
    
    % save sounds
    savename_cur = [savename,'_D_',num2str(D),'_L_',num2str(l),'_'];
    yscale = max(abs([y;ynew;yvoc]));
    wavwrite(0.95*y/yscale,fs,[savedir,'sounds/',savename_cur,'original.wav'])
    wavwrite(0.95*ynew/yscale,fs,[savedir,'sounds/',savename_cur,'chimera1.wav'])
    wavwrite(0.95*yvoc/yscale,fs,[savedir,'sounds/',savename_cur,'chimera2.wav'])

    save([savedir,savename,'.mat'])
  end
end


figure
hold on
plot(Ds,snr_A,'-k')
plot(Ds,snr_Atrain,'-r')
legend('test error','training error')