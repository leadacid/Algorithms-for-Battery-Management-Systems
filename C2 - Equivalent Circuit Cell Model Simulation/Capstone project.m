
% First, make sure that the ESC toolbox functions are in the path
addpath readonly

% GRADED FUNCTION (do not modify this line)

% function rcValues = tuneModel
%
% This function specifies the resistor and capacitor values that
% you choose to use for either an Rint model or a Thévenin model 
% (see lesson 2.1.3), or an "extended Thévenin model" having 
% additional parallel resistor-capacitor branches.
% 
% If you wish to create an Rint model, simply set rcValues equal
% to the value of R0 (in milliohms).
% If you wish to create a Thévenin model, define rcValues to be 
% a vector having three elements. The first element is R0 (in
% milliohms), the second element is R1 (in milliohms), and the 
% third element is C1 (in kilofarads).
% If you wish to create an extended Thévenin model having 
% additional resistor-capacitor branches, then define rcValues to
% be a vector where the first element is R0 in milliohms, the 
% even-index elements (rcValues(2:2:end)) are resistor values for
% the resistor-capacitor branches, in milliohms, and the remaining
% odd-index elements (rcValues(3:2:end)) are capacitor values for
% the resistor-capacitor branches, in kilofarads.
%
% It is possible to receive full credit for this assignment with
% a well-tuned Thévenin model, but it is also possible to get a 
% much better fit to the data using an extended Thévenin model.

function rcValues = tuneModel

  % BEGIN MODIFYING CODE AFTER THIS
  %rcValues = 21; % Rint Model - 14.3368 mV
  %rcValues = [15.6,13.1,64.7] % Thévenin Model - 4.31771 mV
  %rcValues = [15.55, 13.4, 88.6, 1.70, 74.5] % Thévenin model having two resistor-capacitor branches - 3.66018 mV
  rcValues = [15.45, 12.9, 91.9, 1.70, 80, 0.5, 100] % Thévenin model having three resistor-capacitor branches
                                                     % 3.55139 mV
  rcValues = rcValues(:); % output must be a column vector, not a row vector
end  

% Load data file to use for this capstone project
% - these data have the fields: time, current, voltage
% Also load a skeleton ESC model structure
% - has total capacity and OCV relationship
load readonly/pulseData.mat; % load data to use for this project
load readonly/pulseModel.mat; % load SOC/OCV relationships, capacity
T = 25;                            % Test temperature
deltaT = 1;                        % sampling period for data
Q = getParamESC('QParam',T,model); % total capacity of cell

tk = pulseData.time;    % testing time
ik = pulseData.current; % testing current
vk = pulseData.voltage; % testing voltage

% Get tuning values from user-modified function
rcValues = tuneModel;

% Simulate cell model using user-modified rcValues
R0 = rcValues(1)/1000; % convert milliohms to ohms
R = rcValues(2:2:end)/1000; % convert these also
C = rcValues(3:2:end)*1000; % convert kF to F
RCfact = exp(-deltaT./(R.*C));
    
% Simulate the dynamic states of the model
iRk = zeros(length(RCfact),1); % initial resistor currents
vCk = 0*pulseData.voltage; % initialize capacitor voltages
if ~isempty(RCfact)
  for k = 2:length(vCk)
    iRk = diag(RCfact)*iRk + (1-RCfact)*ik(k-1); % update resistor current
    vCk(k) = R'*iRk; % compute capacitor voltage
  end
end

% Simulate SOC state
z0 = SOCfromOCVtemp(pulseData.voltage(1),25,model); 
zk = z0-cumsum([0;ik(1:end-1)])*deltaT/(Q*3600); 

% Compute voltage estimate
vest = OCVfromSOCtemp(zk,25,model) - vCk - ik.*R0;

% Compare against measured voltage, compute RMSE in mV
rmse = 1000*sqrt(mean((vk - vest).^2));

% Plot some results 
subplot(2,1,1); plot(tk/60,vk,tk/60,vest); 
title('Voltage estimation'); grid on
xlabel('Time (min)'); ylabel('Voltage (V))'); 
legend('True measured voltage','Model voltage','location','southeast');

subplot(2,1,2); plot(tk/60,1000*(vk-vest)); 
title('Voltage estimation errors');
xlabel('Time (min)'); ylabel('Voltage error (mV)');
grid on

% Compute the prospective grade
gradeTable = [15 13 11.5 10 9 8 7 6 5 4.5];
ind = find(rmse<gradeTable,1,'last');
if isempty(ind)
  grade = 0;
else
  grade = ind;
end
fprintf('Your tuning values produced an RMS voltage-prediction error of %g mV\n',rmse);
fprintf('Your grade is calculated to be %d/10\n',grade);
