
% add the toolbox code to Octave's path
addpath readonly
% Now, execute the "processDynamic.m" code through the "runProcessDynamic.m" 
% framework (which loads all necessary data files, executes processDynamic.m
% and saves results
% The code will process data from a P14 cell collected at test temperatures
% of 5, 25, and 45 degC. It will save results to "P14model.mat"
runProcessDynamic
% Note that this code takes some time to execute -- perhaps a few minutes
% You will know that the code has finished running when it displays the
% message "Dynamic model created!"

% Load model data file
load P14model.mat
temps = model.temps % display the temperatures processed to create the model
capacities = model.QParam % display measured total capacities at each temperature
hysteresis = model.GParam

% Plot model fields: resistances
subplot(1,2,1); plot(model.temps,1000*model.R0Param)
xlabel('Temperature (degC)'); ylabel('Resistance (mOhms)'); title('Series resistance versus temperature')
subplot(1,2,2); plot(model.temps,1000*model.RParam)
xlabel('Temperature (degC)'); ylabel('Resistance (mOhms)'); title('R-C resistance versus temperature')
