function figh = plot_disc(ynew,y,fs);

% function figh = plot_disc(ynew,y,fs);
%
% plot the original signal and the discrepancy in the time domain

T = length(y);
ts = [1:T]/fs;

figh = figure;
hold on
plot(ts,y,'-','color',[1,1,1]*0.7); 
plot(ts,y-ynew,'-k');
legend('signal','discrepancy')
xlabel('time /s')
ylabel('signals')
