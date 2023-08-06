
% Set up data for problem... use dataset from lesson 4.2.5 for comparison purposes
load readonly/Lesson425data.mat             % Remember that the true capacity is 10.0 Ah

% First, implement FMWLS
Qnom = 10;                                  % Be as generous as possible by using perfect init
gamma = 0.9;                                % Fading-memory rate constant. Default = 1.0
c1 = 1/Sigmay(1);                           % Correct initialization of c1 recursive value
c2 = Qnom/Sigmay(1);                        % Correct initialization of c2 recursive value
for k = 1:length(x),
  c1 = gamma*c1 + x(k)^2/Sigmay(k);         % Update c1 recursive parameter
  c2 = gamma*c2 + x(k)*y(k)/Sigmay(k);      % Update c2 recursive parameter
end
Qhat = c2/c1;                               % The WLS estimate of capacity Q
Sigmaq = 1/c1;
fprintf('FMWLS estimate of Q (true value = 10.0): %f\n',Qhat);
fprintf('FMWLS 3-sigma bounds on Q: %f\n\n',3*sqrt(Sigmaq));

% Now, implement WTLS -- not recursive, uses all data
sigmax = sqrt(Sigmax); % std-dev of x 
sigmay = sqrt(Sigmay); % std-dev of y
Qhat = 5; % purposefully bad initialization of Qhat
for k = 1:10,
  jacobian = sum((2*(Qhat*x-y).*(Qhat*y.*sigmax.^2+x.*sigmay.^2))./((Qhat^2*sigmax.^2+sigmay.^2).^2));
  hessian = sum((2*sigmay.^4.*x.^2+sigmax.^4.*(6*Qhat^2*y.^2-4*Qhat^3*x.*y) - ...
                 sigmax.^2.*sigmay.^2.*(6*Qhat^2*x.^2-12*Qhat*x.*y+2*y.^2))./((Qhat^2*sigmax.^2+sigmay.^2).^3));
  Qhat = Qhat - jacobian/hessian;
  fprintf('WTLS estimate after %d Newton-Raphson iterations: %f\n',k,Qhat);
end
fprintf('WTLS 3-sigma bounds on Q: %f\n\n',3*sqrt(2/hessian));
