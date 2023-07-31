
% First, generate 5000 completely uncorrelated random 2x1 vectors having Gaussian distribution with zero mean
x = randn(2,5000);

% Next, specify the covariance and mean you desire for the output random variables
Sigma = [2 0; 0 0.5]; % covariance matrix for output random variable
mean = [1; 1];          % mean vector for output random variable

% Generate output random variables "y" having desired mean and covariance
% In MATLAB, we would use [L,D] = ldl(Sigma); A = sqrt(D)*L'; 
% But, Octave does not have the ldl function. Instead, we use the LU decomposition
[L,U]=lu(Sigma); A = diag(sqrt(diag(U)))*L';
x = mean(:,ones([1 5000]))+A'*x;

% Plot the output random variables as a scatter plot to view the results
plot(x(1,:),x(2,:),'.');
axis([-5 5 -5 5]); axis('square'); grid on;
title('Scatter plot'); xlabel('x1'); ylabel('x2');
