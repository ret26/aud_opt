clear; 

% compares the normal squared error-objective to one which uses
% average SNR in each channel as the objective. The average SNR
% measure appears harder to optimise and produces worse results.
%
% even intialising using the results of the standard objective did not
% appear to improve things -- perhaps softening the log might help?
% I.e. use log(squared_error+delta)


randn('state',1);

% load sound
[y,fs] = wavread('/home/rich/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/74 - Sentences.wav');

% pick a short section
% y = y(8438:12190);
% y = y(8438+3887:8640+8438);

% pick a long section
 y = y(8438:33630);

% normalise
y = y/sqrt(var(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set auditory model 
DS = 1; % downsampling (optional)

T=length(y);
% Number of channels, slightly less than 1 ERB(Cambridge) per channel.

%channels_per_erb = 0.2; D=ceil(freqtoerb(fs/2)*channels_per_erb); %fc=erbspace(50,fs/2-1000,D);

%D = 20;   fc=erbspace(200,5500,D); bet = 1;
D = 10;   fc=erbspace(200,5500,D); bet = 1;
%D = 3;   fc=erbspace(200,5500,D); bet = 3;
%D = 4;   fc=erbspace(200,5500,D); bet = 2;
%D = 5;   fc=erbspace(200,5500,D); bet = 2;

%D = 3; fc = [100,500,4000]; bet = 3;
%D = 2; fc = [200,800]; bet = 5;
%D = 1; fc = [100]; bet = 7;

filterlength = 5000; % not sure whether this is useful -- try optimising later
betas = bet*ones(1,D); % optional widening of the filters

% choose gammatone filter coefficients
g_gam=gammatonefir(fc,fs,filterlength,betas,'peakphase');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot filter bank
hFig = plotFB(fc,betas,g_gam,fs,T,DS);
drawnow;
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keyboard

% hair cell
% (soft) half wave rectify
rho = 1e4;
fCutLP = 770*2/fs;
ordLP = 7;

% compute targets 
[ATar,YHWTar,YTar] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimise new signal

yInit = 0*randn(T,1)/100;
numIts = ones(10,1)*80;
%numIts = ones(40,1)*40;

[y1,info1] = match_envelopes_V1(yInit,ATar,g_gam,DS,fCutLP,ordLP,rho,numIts,y);
[y2,info2] = match_envelopes_V1(yInit,ATar,g_gam,DS,fCutLP,ordLP,rho,numIts,y,1,0,'log');


[A1,YHW1,Y1] = aud_mod_V1(y1,g_gam,DS,fCutLP,ordLP,rho);
[A2,YHW2,Y2] = aud_mod_V1(y2,g_gam,DS,fCutLP,ordLP,rho);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute similarities to original signal

[snr_y1,snr_A1] = compute_sig_sim(y,y1,fs);
[snr_y2,snr_A2] = compute_sig_sim(y,y2,fs);

snr_Atrain1 = snr(ATar,A1);
snr_Atrain2 = snr(ATar,A2);

disp(['SNR reconstruction (lin) ',num2str(snr_y1),'dB'])
disp(['SNR reconstruction envelopes (lin) ',num2str(mean(snr_A1)),'dB'])
disp(['SNR reconstruction on training envelopes (lin) ',num2str(mean(snr_Atrain1)),'dB'])

disp(['SNR reconstruction (log) ',num2str(snr_y2),'dB'])
disp(['SNR reconstruction envelopes (log) ',num2str(mean(snr_A2)),'dB'])
disp(['SNR reconstruction on training envelopes (log) ',num2str(mean(snr_Atrain2)),'dB'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure
subplot(2,2,1)
hold on
ylabel('objective, lin')
plot(info1.obj)
set(gca,'yscale','log')

subplot(2,2,2)
hold on
plot(info1.err_y)
ylabel('error y, lin')
set(gca,'yscale','log')

subplot(2,2,3)
hold on
ylabel('objective, log')
plot(info2.obj)
set(gca,'yscale','log')

subplot(2,2,4)
hold on
plot(info2.err_y)
ylabel('error y, log')
set(gca,'yscale','log')


figure
hold on
plot(snr_A1,'-k')
plot(snr_A2,'-r')
legend('lin','log')
xlabel('channel number')
ylabel('SNR envelopes')

figure
subplot(2,1,1)
hold on
disc = y - y1;
spec_disc = abs(fft(disc));
spec_y = abs(fft(y));
spec_y1 = abs(fft(y1));

freq = linspace(0,fs/2,floor(T/2));
plot(freq,spec_y(1:floor(T/2)),'-k')
plot(freq,spec_y1(1:floor(T/2)),'-r')
plot(freq,spec_disc(1:floor(T/2)))
ylabel('spectrum, lin')
xlabel('frequency')

subplot(2,1,2)
hold on
disc = y - y2;
spec_disc = abs(fft(disc));
spec_y2 = abs(fft(y2));

plot(freq,spec_y(1:floor(T/2)),'-k')
plot(freq,spec_y2(1:floor(T/2)),'-r')
plot(freq,spec_disc(1:floor(T/2)))
ylabel('spectrum, log')
xlabel('frequency')

figure
subplot(2,1,1)
hold on
plot(y,'-','color',[1,1,1]*0.7); 
plot(y-y1,'-k');
legend('signal','discrepancy')
xlabel('time /samples')
ylabel('signals, lin')

subplot(2,1,2)
hold on
plot(y,'-','color',[1,1,1]*0.7); 
plot(y-y2,'-k');
legend('signal','discrepancy')
xlabel('time /samples')
ylabel('signals, log')


figure
subplot(2,1,1)
hold on
plot(y,'-k')
plot(y1,'-r')
xlabel('time /samples')
ylabel('signals')
legend('original','optimised (lin)')

subplot(2,1,2)
hold on
plot(y,'-k')
plot(y2,'-r')
xlabel('time /samples')
ylabel('signals')
legend('original','optimised (log)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save sounds
yscale = max(abs([y;y1;y2]));
savefile = '/home/rich/Data/aud_opt/demo_compare_log_obj/';
savebasename = ['sentence_bet_',num2str(bet),'_D_',num2str(D),'_']
wavwrite(0.95*y/yscale,fs,[savefile,savebasename,'original.wav'])
wavwrite(0.95*y1/yscale,fs,[savefile,savebasename,'chimera_lin.wav'])
wavwrite(0.95*y2/yscale,fs,[savefile,savebasename,'chimera_log.wav'])
