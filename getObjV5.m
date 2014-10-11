function [obj,dobj] = getObjV5(y,YHWtar,g_gam,DS,rho)

% % Analysis transform
Y = ufilterbank(y,g_gam,DS);
Y = real(Y);

[T,D] = size(Y);

% hair cell
% (soft) half wave rectify
YHW = softHWR(Y,rho);

% objective
dYHW = YHW - YHWtar;
obj1 = 1/2*sum(dYHW(:).^2);

% derivatives

%dAdYHW = filter(z,p,dA);

% dyHWdy = 1./(exp(-rho*y)+1);
% dobj1 = filter(z,p,dA(end:-1:1));
% dobj1 = dobj1(end:-1:1);
% dobj1 = dobj1.*dyHWdy;


% derivative of the soft rectification wrt to input to soft rect
dYHWdY = 1./(exp(-rho*Y)+1);

dObjdY = dYHW.*dYHWdY;

dYHWdy = zeros(T,D);
for d=1:D
  dYHWdy(:,d) = ufilterbank(dObjdY(end:-1:1,d),g_gam(d),DS);
end

dYHWdy = real(dYHWdy(end:-1:1,:));

obj = obj1/T;
dobj = sum(dYHWdy,2)/T;
