
% Don't forget to < shift >< enter> to execute a code cell!
addpath readonly
runProcessOCV

% The "model" structure has three fields of interest to us here: model.SOC, model.OCV0, and model.OCVrel.
% To compute the OCV relationship at temperature T, you would execute:
T = -5; % for example. Change to whatever temperature is of interest to you
OCV_T = model.OCV0 + T*model.OCVrel;
plot(model.SOC,OCV_T); % for example, to visualize results

% We might want to compute the OCV at some desired SOC.
testSOC = 0.25; % for example. Change "0.85" to the SOC of interest
testOCV = interp1(model.SOC,OCV_T,testSOC)
