function figh = plot_compare_sigs(ynew,y,fs);

% function figh = plot_compare_sigs(ynew,y,fs);
%
% plots signals in the time domain

T = length(y);
ts = [1:T]/fs;

figh=figure;
hold on
plot(ts,y,'-k')
plot(ts,ynew,'-r')
xlabel('time /samples')
ylabel('signals')
legend('original','optimised')
