
% Initialize simulation variables
SigmaW = 1; % Process noise covariance
SigmaV = 2; % Sensor noise covariance
maxIter = 40;

% Seed the random number generator: Octave's "randn" function produces pseudo 
% random numbers having a Gaussian distribution. To get the same random numbers
% every time you run the code, you can "seed" the pseudo random number generator 
% with a deterministic value. This allows us to get reproducible results that 
% still contain apparent randomness.
%
% YOU MUST NOT CHANGE THE FOLLOWING LINE IN ORDER TO ARRIVE AT "CORRECT" ANSWERS 
% FOR THE QUIZ. IF YOU CHANGE IT TO EXPLORE THE BEHAVIOR OF THE FILTER FOR 
% DIFFERENT PSEUDO-RANDOM SEQUENCES, CHANGE IT BACK TO randn("seed",-1) BEFORE 
% COMPUTING RESULTS TO ENTER AS QUIZ RESPONSES.
randn("seed",-1);

% Initialize true state, state estimate, error covariance, initial input
xtrue = 2 + randn(1);  % Initialize true system initial state
xhat = 2;              % Initialize Kalman filter initial estimate
SigmaX = 1;            % Initialize Kalman filter covariance
u = 0;                 % Unknown initial driving input: assume zero

% Reserve storage for variables we might want to plot/evaluate
xstore = zeros(maxIter+1,length(xtrue)); xstore(1,:) = xtrue;
xhatstore = zeros(maxIter,length(xhat));
SigmaXstore = zeros(maxIter,length(xhat)^2);

for k = 1:maxIter,
  % EKF Step 0: Compute Ahat, Bhat
  % Note: For this example, x(k+1) = sqrt(5+x(k)) + w(k)
  Ahat = 0.5/sqrt(5+xhat); Bhat = 1;

  % EKF Step 1: State estimate time update
  % Note: You need to insert your system's f(...) equation here
  xhat = sqrt(5+xhat); 

  % KF Step 2: Error covariance time update
  SigmaX = Ahat*SigmaX*Ahat' + Bhat*SigmaW*Bhat';

  % [Implied operation of system in background, with
  % input signal u, and output signal z]
  w = chol(SigmaW)'*randn(1);
  v = chol(SigmaV)'*randn(1);
  ztrue = xtrue^3 + v;  % z is based on present x and u
  xtrue = sqrt(5+xtrue) + w;  % future x is based on present u

  % KF Step 3: Estimate system output
  % Note: You need to insert your system's h(...) equation here
  Chat = 3*xhat^2; Dhat = 1;
  zhat = xhat^3;

  % KF Step 4: Compute Kalman gain matrix
  L = SigmaX*Chat'/(Chat*SigmaX*Chat' + Dhat*SigmaV*Dhat');

  % KF Step 5: State estimate measurement update
  xhat = xhat + L*(ztrue - zhat);
  xhat = max(-5,xhat); % don't get square root of negative xhat!

  % KF Step 6: Error covariance measurement update
  SigmaX = SigmaX - L*Chat*SigmaX;

  % [Store information for evaluation/plotting purposes]
  xstore(k+1,:) = xtrue;
  xhatstore(k,:) = xhat;
  SigmaXstore(k,:) = SigmaX(:);
end

subplot(1,2,1);
plot(0:maxIter-1,xstore(1:maxIter),'k-',0:maxIter-1,xhatstore,'b--', ...
  0:maxIter-1,xhatstore+3*sqrt(SigmaXstore),'m-.',...
  0:maxIter-1,xhatstore-3*sqrt(SigmaXstore),'m-.'); grid;
legend('true','estimate','bounds');
title('Extended Kalman filter in action');
xlabel('Iteration'); ylabel('State');

subplot(1,2,2)
estErr = xstore(1:maxIter)-xhatstore; 
bounds = 3*sqrt(SigmaXstore);
plot(0:maxIter-1,estErr,'b-',0:maxIter-1, bounds,'m--',0:maxIter-1,-bounds,'m--');
grid; legend('Error','bounds',0);
title('EKF Error with bounds');
xlabel('Iteration'); ylabel('Estimation Error');

% Compute rms error
sqrt(mean(estErr.^2))

% Compute fraction of estimates outside of 3-sigma bounds
length(find(abs(estErr)>bounds))/maxIter

max(abs(estErr))
