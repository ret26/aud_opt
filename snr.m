function x = snr(yClean,yNoisy)

% function x = snr(yClean,yNoisy)
%
% Computes the signal to noise ratio between a clean and noisy
% vector:
%
% snr = 10*(log10(sum(yClean.^2))-log10(sum((yClean-yNoisy).^2)))
% 
% INPUTS
% yClean = clean signal size [N,1]
% yNoisy = noisy signal size [N,1]
%
% OUTPUTS
% snr = signal to noise ratio

x = 10*(log10(sum(yClean.^2))-log10(sum((yClean-yNoisy).^2)));