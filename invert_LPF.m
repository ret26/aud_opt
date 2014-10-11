function YHWrecon = invert_LPF(A,rho,fCutLP,ordLP,max_gain)

% function YHWrecon = invert_LPF(A,rho,fCutLP,ordLP,max_gain)
%
% inverts the low pass filter at the end of the model

[T,D] = size(A);
[z,p] = butter(ordLP,fCutLP);

% test to invert filter
yDelta = zeros(2*T,1); % used the wrapping trick to avoid artifacts
yDelta(1) = 1;
aImp = filter(z,p,yDelta);

% get the FFT of the butterworth filter
spec = fft(aImp);
spec_old = abs(spec); 

% invert the fft based filter
spec_new = 1./abs(spec);
spec_new(spec_new>max_gain)= max_gain;

YHWrecon = zeros(T,D);
for d=1:D
  filt_chan = ifft(fft([A(:,d);A(end:-1:1,d)]).*spec_new.*exp(-i*angle(spec)));
  filt_chan = real(filt_chan(1:T));
  YHWrecon(:,d) = filt_chan;
  ind = filt_chan<0;
  YHWrecon(ind,d) = 0;
end


