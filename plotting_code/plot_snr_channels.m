function figh = plot_snr_channels(ynew,y,fc,fs,g_gam,DS,fCutLP,ordLP,rho);

% training envelopes
[A1,YHW,Y] = aud_mod_V1(ynew,g_gam,DS,fCutLP,ordLP,rho);
[A1Tar,YHW,Y] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho);


% test envelopes
channels_per_erb = 1;
filterlength = 5000; % not sure whether this is useful -- try optimising later

T=length(y);
D=ceil(freqtoerb(fs/2)*channels_per_erb);
fc_test=erbspace(50,fs/2-1000,D);
betas = 1*ones(1,D); 
g_gam_test=gammatonefir(fc_test,fs,filterlength,betas,'peakphase');

rho_test = 1e4;
fCutLP_test = 770*2/fs;
ordLP_test = 7;

[A2,YHW2Tar,Y2Tar] = aud_mod_V1(ynew,g_gam_test,1,fCutLP_test, ...
			      ordLP_test,rho_test);
[A2Tar,YHW2,Y2] = aud_mod_V1(y,g_gam_test,1,fCutLP_test, ...
				 ordLP_test,rho_test);


snr_A1 = snr(A1Tar,A1);
snr_A2 = snr(A2Tar,A2);



figh = figure;
hold on
plot(fc,snr_A1,'-k','linewidth',2)
plot(fc_test,snr_A2,'-r','linewidth',2)
xlabel('channel frequency')
ylabel('SNR envelopes')
legend('training SNR','testing SNR')