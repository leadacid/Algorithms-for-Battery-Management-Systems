
% load the E2 circuit model
addpath readonly
load readonly/E2model.mat

% Plot the OCV relationship at a specific temperature
soc = 0:0.01:1; % make a vector of state-of-charge points between 0 and 1 in 1% increments
ocv = OCVfromSOCtemp(soc,25,model);
plot(100*soc,ocv); % "100*" converts to percent
xlabel('State of charge (%)'); ylabel('Voltage (V)'); title('OCV versus SOC relationship for E2 cell')

% Notice how OCVfromSOCtemp.m extrapolates the relationship linearly for out-of-bound states of charge
soc2 = [-0.05:0.01:0, NaN, 1:0.01:1.1]; % NaN to put a break in the plot
ocv2 = OCVfromSOCtemp(soc2,25,model);
plot(100*soc,ocv,'b--',100*soc2,ocv2,'r');
xlabel('State of charge (%)'); ylabel('Voltage (V)'); title('OCV versus SOC relationship for E2 cell')
legend('Normal SOC range','Out-of-bounds SOC range','location','southeast')

% Plot the per-degree change in OCV 
soc = 0:0.01:1;
OCV0 = OCVfromSOCtemp(soc,0,model); % OCV at 0 degrees C
OCV1 = OCVfromSOCtemp(soc,1,model); % OCV at 1 degree C
OCVdiff = 1000*(OCV1-OCV0); % 1000* converts V to mV
plot(100*soc,OCVdiff); % 100* converts to percent
xlabel('State of charge (%)'); ylabel('Voltage change (mV)'); title('Per-degree OCV change for E2 cell');

% Plot the inverse relationship at a specific temperature
ocv = 2:0.01:4.3;
soc = SOCfromOCVtemp(ocv,25,model);
plot(ocv,100*soc); % 100* converts to percent
ylabel('State of charge (%)'); xlabel('Voltage (V)'); title('SOC versus OCV relationship for E2 cell')

ocv = OCVfromSOCtemp(0.25,25,model);


ocv

OCV0 = OCVfromSOCtemp(0.25,0,model); % OCV at 0 degrees C
OCV1 = OCVfromSOCtemp(0.25,1,model); % OCV at 1 degree C
OCVdiff = 1000000*(OCV1-OCV0); % 1000* converts V to mV

OCVdiff
