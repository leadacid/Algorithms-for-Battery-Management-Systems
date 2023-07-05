
% load the E2 circuit model 
addpath readonly
load readonly/E2model.mat

% To retreive a parameter value at a specific temperature, use syntax 
% like this:
getParamESC('R0Param',40,model)


% You can also plot results:
T = 0:1:40;
R0data = zeros(size(T));
for k = 1:length(T),
  R0data(k) = getParamESC('R0Param',T(k),model);
end
plot(T,R0data);

R0data
