
%% Set up data for problem
load readonly/Lesson425data.mat

%% FMTLS
SY = Sigmay(1);  % Assume always this value (which is true in this dataset)
SX = Sigmax(1);  % Assume always this value (which is true in this dataset)
K = sqrt(SX/SY); % Constant of proportionality

Qnom = 10;                        % Initialize nominal capacity
gamma = 1;                     % Fading-memory rate constant. Default = 1.00
SYinit = 1e-3;                    % Uncertainty of Q with respect to Qnom
c1 = 1/SYinit;                    % Correct initialization of c1 recursive value
c2 = Qnom/SYinit;                 % Correct initialization of c2 recursive value
c3 = Qnom^2/SYinit;               % Correct initialization of c3 recursive value

Qhat = 0*x;                       % Initialize storage for output estimate
SigmaQ = 0*x;                     % Initialize storage for estimation-error variance
for k = 1:length(x),
  c1 = gamma*c1 + x(k)^2/SY;      % Update c1 recursive parameter
  c2 = gamma*c2 + x(k)*y(k)/SY;   % Update c2 recursive parameter
  c3 = gamma*c3 + y(k)^2/SY;      % Upcate c3 recursive parameter
  
  Qhat(k) = (-c1+K^2*c3+sqrt((c1-K^2*c3)^2+4*K^2*c2^2))/(2*K^2*c2);
  H = ((-4*K^4*c2)*Q^3+(6*K^4*c3-6*c1*K^2)*Q^2+12*c2*K^2*Q+2*(c1-K^2*c3))/(Q^2*K^2+1)^3;
  
  SigmaQ(k) = 2/H;
end

xvals = 1:length(Qhat);
plot(xvals,10*ones(size(Qhat)),xvals,Qhat,[xvals,NaN,xvals],[Qhat+3*sqrt(SigmaQ),NaN,Qhat-3*sqrt(SigmaQ)]);
xlabel('Update iteration number'); ylabel('Total capacity estimate (Ah)');
title('Output of FMTLS method');   xlim([0 1000])
legend('True total capacity','FMTLS estimate','Bounds on estimate','location','northeast');

Qhat

3*sqrt(SigmaQ)
