function [A,YHW,Y] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho)

% function [A,YHW,Y] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho)
%
% implements simple auditory model with gammatone filterbank,
% threshold and low pass filter auditory nerve model
%
% INPUTS
% y = signal [T,1]
% g_gam = gammatone filter coefficients produced from gammatonefir []
% DS = down sampling option, not currently implemented
% fCutLP = low pass filter cutoff for AN model
% ordLP = order of the butterworth implementation
% rho = strength of the soft-threshold-linear function, log(1+exp(rho*x))/rho
%
% OUTPUTS
% A = auditory nerve fibre outputs, [T,D]
% YHW = halve wave rectified gammatone filter outputs [T,D]
% Y = gammatone filter outputs [T,D]

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

