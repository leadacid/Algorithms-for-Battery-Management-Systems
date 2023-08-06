
% First load the dataset: contains vectors of (x,y) data in "x" and "y", and variances of
% the data in "Sigmax" and "Sigmay". Note that the true capacity used to create this dataset was Q=10.0.
load readonly/Lesson425data.mat

% First, implement the WLS method for comparison purposes. 
c1 = sum(x.^2./Sigmay);
c2 = sum(x.*y./Sigmay);
fprintf('WLS estimate of Q (true value = 10.0): %f\n\n',c2/c1);

% Now, implement the WTLS method. 
sigmax = sqrt(Sigmax); % convert variance into standard deviation
sigmay = sqrt(Sigmay); % convert variance into standard deviation
Qhat = 5;              % purposefully bad initialization of Qhat to show that method converges
for k = 1:10,          % Newton-Raphson iterations to update estimate of Qhat
  jacobian = sum((2*(Qhat*x-y).*(Qhat*y.*sigmax.^2+x.*sigmay.^2))./((Qhat^2*sigmax.^2+sigmay.^2).^2));
  hessian = sum((2*sigmay.^4.*x.^2+sigmax.^4.*(6*Qhat^2*y.^2-4*Qhat^3*x.*y) - ...
                 sigmax.^2.*sigmay.^2.*(6*Qhat^2*x.^2-12*Qhat*x.*y+2*y.^2))./((Qhat^2*sigmax.^2+sigmay.^2).^3));
  Qhat = Qhat - jacobian/hessian;
  fprintf('WTLS estimate of Q after %2d Newton-Raphson iterations: %f\n',k,Qhat);
end
