
addpath ./readonly             % Add helper functions to Octave's path
load ./readonly/CellModel.mat  % Load the ESC cell model we will use

% Simulate a cell for K timesteps, with input ik and starting state x0
% x0 = [z0; ir0; h0]
function [vDT,xDT] = simCellKDT(ik,x0,A,B,KDT,T,model,R0,R,M,M0)

  % Form the "A" and "B" matrices needed
  Amat = A(ik); Bmat = B(ik); dA = diag(Amat);
  
  if ik == 0, % AK = matrix to multiply B and input (see lesson)
    ADT = diag([KDT, (1-dA(2)^KDT)/(1-dA(2)), KDT]);
  else
    ADT = diag([KDT, (1-dA(2)^KDT)/(1-dA(2)), (1-dA(3)^KDT)/(1-dA(3))]);
  end
  xDT = (dA).^KDT.*x0 + ADT*Bmat*[ik; sign(ik)]; % efficiently compute future state
  
  % Compute voltage based on future state
  vDT = OCVfromSOCtemp(xDT(1),T,model) + M*xDT(3) + M0*sign(ik) - R*xDT(2) - ik*R0;
end

% Set operating temperature
T = 25;

% Get ESC model parameter values
R = getParamESC('RParam',T,model);
RC = getParamESC('RCParam',T,model);
M = getParamESC('MParam',T,model);
M0 = getParamESC('M0Param',T,model);
R0 = getParamESC('R0Param',T,model);
Gamma = getParamESC('GParam',T,model);
Q = getParamESC('QParam',T,model);

A = @(ik) diag([1 exp(-1/(RC)) exp(-abs(ik*Gamma/(3600*Q)))]);
B = @(ik) [-1/(3600*Q) 0; (1-exp(-1/RC)) 0; ...
           0 (exp(-abs(ik*Gamma/(3600*Q)))-1)];
           
% Initial state set to 50% SOC, at rest (resistor currents=0), no dynamic hysteresis
x0 = [0.5; 0; 0];
[vDT,xDT] = simCellKDT(1,x0,A,B,10,T,model,R0,R,M,M0) % simulate 1A discharge for 10 seconds at 25 degC.
