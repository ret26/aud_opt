function ynew = cosine_ramp(y,Tsmooth);

% function ynew = cosine_ramp(y,Tsmooth);
%
% add cosine smoothing to the start and end of a signal

T = length(y);
sinOn = sin(2*pi*[0:Tsmooth]'/(4*Tsmooth));
cosine_smooth = [sinOn;ones(T-2*Tsmooth-2,1);sinOn(end:-1:1)];
ynew = y.*cosine_smooth;
