function [ynew,info] = match_HWR_filt_V1(ynew,YHWTar,g_gam,DS,rho,numIts,y,varargin)

% function [ynew,info] = match_HWR_filt_V1(ynew,YHWTar,g_gam,DS,rho,numIts,y,varargin)
%
%
% INPUTS
% ynew = initialisation for new signal [T,1]
% YHWTar = target half wave rectified filter coefficient [T,D]
% g_gam = gammatone filter coefficients produced from gammatonefir []
% DS = down sampling option, not currently implemented
% rho = strength of the soft-threshold-linear function, log(1+exp(rho*x))/rho
% numIts = number of iterations [L,1]
% y = signal approximating (only required to give back information
%     about ynew's proximity to the original signal)
% optional arguments:
% verbose = flag to indicate whether to display information about
%    the optimisation

% OUTPUTS
% ynew = optimised signal
% info = structure of algorithmic information including
%      obj = objective function evaluations
%      err_y = mean squared error to signal y

I = length(numIts);
obj = [];
it = [];
err_y = [];

if nargin>9
  verbose = varargin{1};
else
  verbose = 1;
end


for it = 1:I

  if verbose==1
    disp([num2str(it),'/',num2str(I)])
  end
  
  itSet.length = numIts(it);
  itSet.method = 'CG';
%  itSet.mem = 500;

tic;

    [ynew, objCur, itCur] = minimize_new(ynew,'getObjV5',itSet, ...
					 YHWTar,g_gam,DS,rho); 
    timCur = toc;
    
  obj = [obj;objCur];
  it = [it;itCur];

  
  err_y_cur = mean((ynew-y).^2);
  err_y = [err_y,err_y_cur];
  
  if verbose==1
    disp(['objective ',num2str(objCur(end)),'   time ',num2str(timCur),'s'])
    disp(['error y  ',num2str(err_y(end)),'   time ',num2str(timCur),'s'])
  end
  
end

info.err_y = err_y;
info.obj = obj;
