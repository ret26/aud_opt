function [obj,dobj] = getObjV3(y,Atar,g_gam,DS,fCutLP,ordLP,rho,kappa)

% 
% filter bank: gammatone
% envelopes: half wave rectify
% likelihood cost: squared error
% prior: squared error penalty on y

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

% yHW = softHWR(y,rho);
% [z,p] = butter(ordLP,fCutLP);
% A = filter(z,p,yHW);

% objective
dA = A - Atar;

sum_dAsq = sum(dA.^2,1);

obj1 = sum(log(sum_dAsq/T)); %1/2*sum(dA(:).^2);

obj2 = 1/2*kappa*sum(y.^2);

% derivatives

%dAdYHW = filter(z,p,dA);

% dyHWdy = 1./(exp(-rho*y)+1);
% dobj1 = filter(z,p,dA(end:-1:1));
% dobj1 = dobj1(end:-1:1);
% dobj1 = dobj1.*dyHWdy;


% derivative of hair-cell low pass output wrt halve w
dAdYHW = zeros(T,D);
for d=1:D
  dAdYHW(:,d) = filter(z,p,dA(end:-1:1,d));
end

dAdYHW = dAdYHW(end:-1:1,:);

% derivative of the soft rectification wrt to input to soft rect
dYHWdY = 1./(exp(-rho*Y)+1);

dAdY = dAdYHW.*dYHWdY;

dAdy = zeros(T,D);
for d=1:D
  dAdy(:,d) = ufilterbank(dAdY(end:-1:1,d),g_gam(d),DS);
end

dAdy = real(dAdy(end:-1:1,:));

%putting this together
obj = obj1 + obj2/T;
dobj = sum(dAdy./(ones(T,1)*sum_dAsq),2)*2 + kappa*y/T;


