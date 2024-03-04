% --------------------------------------------------------------------
% script runProcessDynamic
%
% RUNPROCESSDYNAMIC reads data files corresponding to dynamic cell 
% tests, executes PROCESSDYNAMIC, and then saves the resulting 
% model.  It relies on SETUPDYNDATA to provide a list of data files
% to be processed.

setupDynData; % get list of files to be processed
numpoles = 1; % number of resistor--capacitor pairs in final model
doHyst = 1;   % whether to include hysteresis in model

for indID = 1:length(cellIDs), % process each cell type
  cellID = cellIDs{indID};
  
  % Read model OCV file
  modelFile = sprintf('%smodel-ocv.mat',cellID);
  load(modelFile);
  
  % Read MAT raw data files
  data = zeros([0 length(mags{indID} > 0)]); dataInd = 0;
  for indTemps = 1:length(mags{indID}),
    theMag = mags{indID}(indTemps);
    if theMag < 0, 
      continue 
    else
      dataInd = dataInd + 1;
    end
    if temps(indTemps) < 0,
      DYNPrefix = sprintf('readonly/%s_DYN/%s_DYN_%02d_N%02d',...
        cellID,cellID,theMag,abs(temps(indTemps)));
    else        
      DYNPrefix = sprintf('readonly/%s_DYN/%s_DYN_%02d_P%02d',...
        cellID,cellID,theMag,temps(indTemps));
    end
    inFile = sprintf('%s.mat',DYNPrefix);
    fprintf('Loading %s\n',inFile);
    load(inFile);        
    data(dataInd).temp    = temps(indTemps); % temperature
    data(dataInd).script1 = DYNData.script1;
    data(dataInd).script2 = DYNData.script2;
    data(dataInd).script3 = DYNData.script3;
  end
  
  model = processDynamic(data,model,numpoles,doHyst);
  modelFile = sprintf('%smodel.mat',cellID);
  save(modelFile,'model');

  fprintf('\nDynamic model created!\n');
end