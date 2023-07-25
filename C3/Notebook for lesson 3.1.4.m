
% Add model and data files, plus ESC toolbox to path
addpath ./readonly
load('readonly/CellData.mat');
load('readonly/CellModel.mat');

socEstimate1 = SOCfromOCVtemp(voltage,25,model);

% SOC Estimation Plot
t = 0:length(socEstimate1)-1; t = t/60;
plot(t,100*socEstimate1, t,100*soc);
title('True SOC and voltage-based estimate');
xlabel('Time (min)'); ylabel('SOC and estimate (%)');
axis([0 500 0 100]); legend('SOC estimate','True SOC'); grid on;

% Compute answer for quiz based on first kind of SOC estimation here...
sqrt(mean((soc-socEstimate1).^2)) * 100

R0=getParamESC('R0param',25,model);

socEstimate2 = SOCfromOCVtemp(voltage+current*R0,25,model);

% SOC Estimation Plot
t = 0:length(socEstimate2)-1; t = t/60;
plot(t,100*socEstimate2, t,100*soc);
title('True SOC and voltage-based estimate');
xlabel('Time (min)'); ylabel('SOC and estimate (%)');
axis([0 500 0 100]); legend('SOC estimate','True SOC'); grid on;

% Compute answer for quiz based on second kind of SOC estimation here...
sqrt(mean((soc-socEstimate2).^2)) * 100
