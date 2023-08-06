
% Set up data for problem: Note that the "x" vector was chosen deterministically for this problem
% and the "y" vector was generated randomly assuming a true capacity of 10 Ah.
Sigmay = 0.2^2;                             % Variance of random noise on "y" when generating data
x = [0 0.2 0.4 0.6 0.8 1.0];                % No noise on "x" in this example
y = [-0.04, 2.18, 3.85, 5.72, 7.72, 10.10]; % y = Q*x + sqrt(Sigmay)*randn(size(x));

c1 = sum(x.^2/Sigmay)                       % The denominator of the WLS computation
c2 = sum(x.*y/Sigmay)                       % The numerator of the WLS computation
Qhat = c2/c1                                % The WLS estimate of capacity Q
