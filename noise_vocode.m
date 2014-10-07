function yvoc2 = noise_vocode(ATar,g_gam,DS,fCutLP,ordLP,rho)

% function yvoc2 = noise_vocode(ATar,g_gam,DS,fCutLP,ordLP,rho)
%
% noise vocode according to Shannon 1995

[T,D] = size(ATar);

yvoc2 = zeros(T,1);
for d=1:D
  yCur = randn(T,1).*ATar(:,d);
  yCur = ufilterbank(yCur,g_gam(d),DS);
  yCur = real(yCur);
  yvoc2 = yvoc2+yCur; 
end

yvoc2 = yvoc2/sqrt(var(yvoc2));
[Avoc2,YHWvoc2,Yvoc2] = aud_mod_V1(yvoc2,g_gam,DS,fCutLP,ordLP,rho);
