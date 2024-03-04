% simCell: Simulate the ESC cell model
%
% [vk,irck,hk,zk,OCV] = simCell(ik,T,model,deltaT,z0,h0,iR0)
%
% Inputs:  ik    = vector input current (A)
%          T     = temperature (degC)
%          model = ESC model structure
%          deltaT= sampling period (s)
%          z0    = initial cell SOC (0..1)
%          h0    = initial hysteresis state (0..1)
%          iR0   = initial RC-filter current state (A)
%
% Outputs: vk    = cell voltage (V)
%          irck  = RC-filter resistor current(s) (A)
%          hk    = hysteresis state (0..1)
%          zk    = cell state-of-charge (0..1)
%          OCV   = cell OCV value (V)
%
% Important note: time vector corresponding to ik assumed to be
%          tk    = (0:length(ik)-1)*deltaT, 
% where dt = sampling period of model.  Same time vector corresponds 
% to output values, so first entry in all outputs corresponds to 
% tk = 0 (state has not yet advanced), next entry corresponds to 
% tk = deltaT; etc.

function [vk,irck,hk,zk,OCV] = simCell(ik,T,model,deltaT,z0,h0,iR0)
  % Force data to be column vector(s) in case user entered data incorrectly
  ik = ik(:); iR0 = iR0(:); 
  N = length(ik); Nr = length(iR0);
  % initialize some outputs
  irck = zeros(N,Nr); hk = zeros(N,1); zk = zeros(N,1); 
  sik = zeros(N,1); OCV = zeros(N,1);
  irck(1,:) = iR0'; hk(1) = h0; zk(1) = z0; sik(1) = 0;
  
  RCfact = exp(-deltaT./abs(getParamESC('RCParam',T,model)))';
  if length(RCfact) ~= Nr,
    error('iR0 does not have the correct number of entries');
  end
  G = getParamESC('GParam',T,model);
  Q = getParamESC('QParam',T,model);
  M = getParamESC('MParam',T,model);
  M0 = getParamESC('M0Param',T,model);
  RParam = getParamESC('RParam',T,model);
  R0Param = getParamESC('R0Param',T,model);
  etaParam = getParamESC('etaParam',T,model);
  
  etaik = ik; etaik(ik<0) = etaParam*ik(ik<0);

  % Simulate the dynamic states of the model
  for k = 2:length(ik),
    irck(k,:) = irck(k-1,:)*diag(RCfact) + (1-RCfact')*etaik(k-1);
  end
  zk = z0-cumsum([0;etaik(1:end-1)])*deltaT/(Q*3600);
  if any(zk>1.1),
    warning('Current may have wrong sign as SOC > 110%');
  end
  
  % Hysteresis stuff
  fac=exp(-abs(G*etaik*deltaT/(3600*Q)));
  for k=2:length(ik),
    hk(k)=fac(k-1)*hk(k-1)-(1-fac(k-1))*sign(ik(k-1));
    sik(k) = sign(ik(k));
    if abs(ik(k))<Q/100, sik(k) = sik(k-1); end
  end
    
  % Compute output equation
  OCV = OCVfromSOCtemp(zk,T,model);
  
  vk = OCV - irck*RParam' - R0Param*ik + M*hk + M0*sik;
end