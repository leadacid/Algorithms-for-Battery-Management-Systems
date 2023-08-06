
% Set up data for problem: Same dataset as used last week
load readonly/Lesson425data.mat

% TLS method implementation
SY = Sigmay(1);  % Assume always this value (which is true in this dataset)
SX = Sigmax(1);  % Assume always this value (which is true in this dataset)
K = sqrt(SX/SY); % Constant of proportionality 

% Compute summation terms
c1 = cumsum(x.^2/SY); % Value of c1 for every iteration number
c2 = cumsum(x.*y/SY); % Value of c2 for every iteration number
c3 = cumsum(y.^2/SY); % Value of c3 for every iteration number

% Compute solution for every iteration in one statement (slide 5)
Qhat = (-c1+K^2*c3+sqrt((c1-K^2*c3).^2+4*K^2*c2.^2))./(2*K^2*c2);

% Plot some results 
xvals = 1:length(Qhat);
plot(xvals,10*ones(size(Qhat)),xvals,Qhat);
xlabel('Update iteration number'); ylabel('Total capacity estimate (Ah)');
title('Output of TLS method');     xlim([0 1000])
legend('True total capacity','TLS estimate','location','northeast');

Qhat
