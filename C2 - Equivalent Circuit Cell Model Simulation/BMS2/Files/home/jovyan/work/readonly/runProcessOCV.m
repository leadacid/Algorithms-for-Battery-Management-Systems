clear all
cellIDs = {'P14'};
% data files for each cell available at these temperatures
temps = {[-25 -15 -5 5 15 25 35 45]};    % P14
% minimum and maximum voltages for each cell, used for plotting results       
minV = [2.50];
maxV = [4.25];

% --------------------------------------------------------------------
% Load raw data from cell tests, then process
% --------------------------------------------------------------------
for theID = 1:length(cellIDs), 
  dirname = cellIDs{theID}; cellID = dirname;
  ind = find(dirname == '_');
  if ~isempty(ind), dirname = dirname(1:ind-1); end
  OCVDir = sprintf('readonly/%s_OCV',dirname);  
  
  filetemps = temps{theID}(:); 
  numtemps = length(filetemps); 
  data = zeros([0 numtemps]);

  for k = 1:numtemps,
    if filetemps(k) < 0,
      filename = sprintf('%s/%s_OCV_N%02d.mat',...
        OCVDir,cellID,abs(filetemps(k)));
    else
      filename = sprintf('%s/%s_OCV_P%02d.mat',...
        OCVDir,cellID,filetemps(k));
    end
    fprintf('Loading OCV data collected for test temperature %d degrees C\n',filetemps(k));
    load(filename);
    data(k).temp = filetemps(k);       
    data(k).script1 = OCVData.script1;
    data(k).script2 = OCVData.script2;
    data(k).script3 = OCVData.script3;
    data(k).script4 = OCVData.script4;
  end

  fprintf('Processing data to create OCV relationship\n');
  model = processOCV(data,cellID,minV(theID),maxV(theID),0);
  save(sprintf('%smodel-ocv.mat',cellID),'model');
  fprintf('Processing complete. Results are available in variable "model"\n');
  fprintf('and are saved in %s\n',sprintf('%smodel-ocv.mat',cellID));
end