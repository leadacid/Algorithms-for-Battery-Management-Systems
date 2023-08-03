
addpath readonly
load readonly/E2model

% You may have forgotten the syntax for retreiving ESC model parameter values. 
% Here is a refresher:
help getParamESC

% As an example, we can retreive the dynamic hysteresis "gamma" value at 25 degC
hystGamma = getParamESC('GParam',25,model)

% You may have forgotten the syntax for retreiving open circuit voltage at a
% particular state of charge and temperature. Here is a refresher:
help OCVfromSOCtemp

% As an example, we can retreive the OCV at 50% SOC and 25 degrees celsius
ocv = OCVfromSOCtemp(0.5,25,model)

Q = getParamESC('QParam',25,model) * 3600

ahk = exp(-2*hystGamma/Q)

RC = getParamESC('RCParam',25,model)

Ark = exp(-1/RC)

SOC = 0:0.01:1;
OCV = OCVfromSOCtemp(SOC,25,model);

z = 0.2

dZ = SOC(2) - SOC(1);  % Find spacing of SOC vector
dUdZ = diff(OCV)/dZ; % Scaled forward finite difference
dOCV = ([dUdZ(1) dUdZ] + [dUdZ dUdZ(end)])/2; % Avg of fwd / bkwd diffs
dOCVz = interp1(SOC ,dOCV ,z); 

dOCVz
