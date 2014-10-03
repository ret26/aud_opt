function spec = gammatone_spectrum(f0,ban,freqs)

% Returns gammatone spectrum with centre-frequencies f0 and
% bandwidths ban at frequencies freqs
%
% INPUTS
% f0 = centre-frequencies [D,1]
% ban = bandwidth
% freqs = frequencies at which to evaluate [N,1]
%
% OUTPUTS
% spec = gammatone spectrum [N,D]


n=4;
a = pi*factorial(2*n-2)*2^(-(2*n-2))/(factorial(n-1)^2);
b = ban/a;
N= length(freqs);
D = length(f0);
spec = zeros(N,D);
  
for d=1:D
  spec(:,d) = abs((1+i*(freqs-f0(d))/b(d)).^(-4) +  ...
		  (1+i*(freqs+f0(d))/b(d)).^(-4));
end
