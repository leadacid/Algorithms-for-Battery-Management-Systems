% --------------------------------------------------------------------
% function processDynamic
%
% Technical note: PROCESSDYNAMIC assumes that specific Arbin test 
% scripts have been executed to generate the input files. 
% "makeMATfiles.m" converts the raw Excel data files into "MAT" format
% where the MAT files have fields for time, step, current, voltage, 
% chgAh, and disAh for each script run.
%
% The results from three scripts are required at every temperature.  
% The steps in each script file are assumed to be:
%   Script 1 (thermal chamber set to test temperature):
%     Step 1: Rest @ 100% SOC to acclimatize to test temperature
%     Step 2: Discharge @ 1C to reach ca. 90% SOC
%     Step 3: Repeatedly execute dynamic profiles (and possibly
%             intermediate rests) until SOC is around 10%
%   Script 2 (thermal chamber set to 25 degC):
%     Step 1: Rest ca. 10% SOC to acclimatize to 25 degC
%     Step 2: Discharge to min voltage (ca. C/3)
%     Step 3: Rest
%     Step 4: Constant voltage at vmin until current small (ca. C/30)
%     Steps 5-7: Dither around vmin
%     Step 8: Rest
%   Script 3 (thermal chamber set to 25 degC):
%     Step 2: Charge @ 1C to max voltage
%     Step 3: Rest
%     Step 4: Constant voltage at vmax until current small (ca. C/30)
%     Steps 5-7: Dither around vmax
%     Step 8: Rest
% 
% All other steps (if present) are ignored by PROCESSDYNAMIC. The time 
% step between data samples must be uniform -- we assume a 1s sample
% period in this code
%
% The inputs:
% - data: An array, with one entry per temperature to be processed. 
%         One of the array entries must be at 25 degC. The fields of 
%         "data" are: temp (the test temperature), script1, 
%         script 2, and script 3, where the latter comprise data 
%         collected from each script.  The sub-fields of these script 
%         structures that are used by PROCESSDYNAMIC are the vectors: 
%         current, voltage, chgAh, and disAh
% - model: The output from processOCV, comprising the OCV model
% - numpoles: The number of R-C pairs in the model
% - doHyst: 0 if no hysteresis model desired; 1 if hysteresis desired
%
% The output:
% - model: A modified model, which now contains the dynamic fields
%         filled in.
function model = processDynamic(data,model,numpoles,doHyst)
  global bestcost
  
  % used by fminbnd later on
  options=optimset('TolX',1e-8,'TolFun',1e-8,'MaxFunEval',100000, ...
    'MaxIter',1e6,'Jacobian','Off'); % for later optimization
  options=optimset('TolX',0.1,'TolFun',1e-2,'MaxFunEval',40, ...
    'MaxIter',20,'Jacobian','Off'); % for later optimization
    
  % ------------------------------------------------------------------
  % Step 1: Compute capacity and coulombic efficiency for every test
  % ------------------------------------------------------------------
  alltemps = [data(:).temp];
  alletas  = 0*alltemps;
  allQs    = 0*alltemps;
  
  ind25 = find(alltemps == 25); 
  if isempty(ind25),
    error('Must have a test at 25degC');
  end
  not25 = find(alltemps ~= 25);
  
  for k = ind25,    
    totDisAh = data(k).script1.disAh(end) + ...
               data(k).script2.disAh(end) + ...
               data(k).script3.disAh(end);
    totChgAh = data(k).script1.chgAh(end) + ...
               data(k).script2.chgAh(end) + ...
               data(k).script3.chgAh(end);
    eta25 = totDisAh/totChgAh; 
    data(k).eta = eta25; alletas(k) = eta25;
    data(k).script1.chgAh = data(k).script1.chgAh*eta25;
    data(k).script2.chgAh = data(k).script2.chgAh*eta25;
    data(k).script3.chgAh = data(k).script3.chgAh*eta25;    

    Q25 = data(k).script1.disAh(end) + data(k).script2.disAh(end) -...
          data(k).script1.chgAh(end) - data(k).script2.chgAh(end);
    data(k).Q = Q25; allQs(k) = Q25;
  end
  eta25 = mean(alletas(ind25));
  
  for k = not25,    
    data(k).script2.chgAh = data(k).script2.chgAh*eta25;
    data(k).script3.chgAh = data(k).script3.chgAh*eta25;
    eta = (data(k).script1.disAh(end) + data(k).script2.disAh(end)+...
           data(k).script3.disAh(end) - data(k).script2.chgAh(end)-...
           data(k).script3.chgAh(end))/data(k).script1.chgAh(end);
    data(k).script1.chgAh = eta*data(k).script1.chgAh;
    data(k).eta = eta; alletas(k) = eta;
    
    Q = data(k).script1.disAh(end) + data(k).script2.disAh(end) - ...
          data(k).script1.chgAh(end) - data(k).script2.chgAh(end);
    data(k).Q = Q; allQs(k) = Q;
  end
  
  model.temps    = unique(alltemps); numTemps = length(model.temps);
  model.etaParam = NaN(1,numTemps);
  model.QParam   = NaN(1,numTemps);
  for k = 1:numTemps,
    model.etaParam(k) = mean(alletas(alltemps == model.temps(k)));
    model.QParam(k)   = mean(allQs(alltemps == model.temps(k)));
  end
  
  % ------------------------------------------------------------------
  % Step 2: Compute OCV for "discharge portion" of test
  % ------------------------------------------------------------------
  for k = 1:length(data),
    etaParam = model.etaParam(k);
    etaik = data(k).script1.current; 
    etaik(etaik<0)= etaParam*etaik(etaik<0);
    data(k).Z = 1 - cumsum([0,etaik(1:end-1)])*1/(data(k).Q*3600); 
    data(k).OCV = OCVfromSOCtemp(data(k).Z(:),alltemps(k),model);
  end
  
  % ------------------------------------------------------------------
  % Step 3: Now, optimize!
  % ------------------------------------------------------------------
  model.GParam  = NaN(1,numTemps); % "gamma" hysteresis parameter
  model.M0Param = NaN(1,numTemps); % "M0" hysteresis parameter
  model.MParam  = NaN(1,numTemps); % "M" hysteresis parameter
  model.R0Param = NaN(1,numTemps); % "R0" ohmic resistance parameter
  model.RCParam = NaN(numTemps,numpoles); % time const.
  model.RParam  = NaN(numTemps,numpoles); % Rk

  for theTemp = 1:numTemps, 
    fprintf('Processing temperature %d\n',model.temps(theTemp));
    bestcost = Inf;
    if doHyst,
      model.GParam(theTemp) = abs(fminbnd(@(x) optfn(x,data,...
                                  model,model.temps(theTemp),...
                                  doHyst),1,250,options));
    else
      model.GParam(theTemp) = 0;
      theGParam = 0;
      optfn(theGParam,data,model,model.temps(theTemp),doHyst);
    end
    [~,model] = minfn(data,model,model.temps(theTemp),doHyst);                          
  end
return

% --------------------------------------------------------------------
% This minfn works for the enhanced self-correcting cell model
% --------------------------------------------------------------------
function cost=optfn(theGParam,data,model,theTemp,doHyst)
  global bestcost 
  
  model.GParam(model.temps == theTemp) = abs(theGParam);
  [cost,model] = minfn(data,model,theTemp,doHyst);
  if cost<bestcost, % update plot of model params for every improvement
    bestcost = cost;
    disp('    The model created for this value of gamma is the best ESC model yet!');
  end
return

% --------------------------------------------------------------------
% Using an assumed value for gamma (already stored in the model), find 
% optimum values for remaining cell parameters, and compute the RMS 
% error between true and predicted cell voltage
% --------------------------------------------------------------------
function [cost,model]=minfn(data,model,theTemp,doHyst)
  alltemps = [data(:).temp];
  ind = find(alltemps == theTemp); numfiles = length(ind);

  xplots = ceil(sqrt(numfiles));
  yplots = ceil(numfiles/xplots);
  rmserr = zeros(1,xplots*yplots);
  
  G = abs(getParamESC('GParam',theTemp,model));
  Q = abs(getParamESC('QParam',theTemp,model));
  eta = abs(getParamESC('etaParam',theTemp,model));
  RC = getParamESC('RCParam',theTemp,model);
  numpoles = length(RC);
  
  for thefile = 1:numfiles;
    ik = data(ind(thefile)).script1.current(:);
    vk = data(ind(thefile)).script1.voltage(:);
    tk = (1:length(vk))-1;
    etaik = ik; etaik(ik<0) = etaik(ik<0)*eta;

    h=0*ik; sik = 0*ik;
    fac=exp(-abs(G*etaik/(3600*Q)));
    for k=2:length(ik),
      h(k)=fac(k-1)*h(k-1)+(fac(k-1)-1)*sign(ik(k-1));
      sik(k) = sign(ik(k));
      if abs(ik(k))<Q/100, sik(k) = sik(k-1); end
    end
    
    % First modeling step: Compute error with model = OCV only
    vest1 = data(ind(thefile)).OCV;
    verr = vk - vest1;
    
    % Second modeling step: Compute time constants in "A" matrix
    np = numpoles; 
    while 1,
      A = SISOsubid(-diff(verr),diff(etaik),np);
      eigA = eig(A); 
      eigA = eigA(eigA == conj(eigA));  % make sure real
      eigA = eigA(eigA > 0 & eigA < 1); % make sure in range
      okpoles = length(eigA); np = np+1;
      if okpoles >= numpoles, break; end
      fprintf('Trying np = %d\n',np);
    end    
    RCfact = sort(eigA); RCfact = RCfact(end-numpoles+1:end);
    RC = -1./log(RCfact);
    % Simulate the R-C filters to find R-C currents
    vrcRaw = zeros(numpoles,length(h));
    for k=2:length(ik),
      vrcRaw(:,k) = diag(RCfact)*vrcRaw(:,k-1) + (1-RCfact)*etaik(k-1);
    end
    vrcRaw = vrcRaw';

    % Third modeling step: Hysteresis parameters
    if doHyst,
      H = [h,sik,-etaik,-vrcRaw]; 
      W = lsqnonneg(H,verr); %  W = H\verr;    
      M = W(1); M0 = W(2); R0 = W(3); Rfact = W(4:end)';
    else
      H = [-etaik,-vrcRaw]; 
      W = H\verr;    
      M=0; M0=0; R0 = W(1); Rfact = W(2:end)';
    end
    ind = find(model.temps == data(ind(thefile)).temp,1);
    model.R0Param(ind) = R0;
    model.M0Param(ind) = M0;
    model.MParam(ind) = M;
    model.RCParam(ind,:) = RC';
    model.RParam(ind,:) = Rfact';
    
    vest2 = vest1 + M*h + M0*sik - R0*etaik - vrcRaw*Rfact';
    verr = vk - vest2;
        
    % Compute RMS error only on data roughly in 5% to 95% SOC
    v1 = OCVfromSOCtemp(0.95,data(ind(thefile)).temp,model);
    v2 = OCVfromSOCtemp(0.05,data(ind(thefile)).temp,model);
    N1 = find(vk<v1,1,'first'); N2 = find(vk<v2,1,'first');
    if isempty(N1), N1=1; end; if isempty(N2), N2=length(verr); end
    rmserr(thefile)=sqrt(mean(verr(N1:N2).^2));    
  end 

  cost=sum(rmserr); 
  fprintf('  RMS error for present value of gamma = %0.2f (mV)\n',cost*1000);
  if isnan(cost), stop, end
return

% A = SISOsubid(y,u,n);
%  Identifies state-space "A" matrix from input-output data.
%     y: vector of measured outputs
%     u: vector of measured inputs 
%     n: number of poles in solution
%           
%     A: discrete-time state-space state-transition matrix.
%                 
%  Theory from "Subspace Identification for Linear Systems
%               Theory - Implementation - Applications" 
%               Peter Van Overschee / Bart De Moor (VODM)
%               Kluwer Academic Publishers, 1996
%               Combined algorithm: Figure 4.8 page 131 (robust)
%               Robust implementation: Figure 6.1 page 169
%
%  Code adapted from "subid.m" in "Subspace Identification for 
%               Linear Systems" toolbox on MATLAB CENTRAL file 
%               exchange, originally by Peter Van Overschee, Dec. 1995
function A = SISOsubid(y,u,n)
  y = y(:); y = y'; ny = length(y); % turn y into row vector
  u = u(:); u = u'; nu = length(u); % turn u into row vector
  i = 2*n; % #rows in Hankel matrices. Typically: i = 2 * (max order)
  twoi = 4*n;           

  if ny ~= nu, error('y and u must be same size'); end
  if ((ny-twoi+1) < twoi); error('Not enough data points'); end

  % Determine the number of columns in the Hankel matrices
  j = ny-twoi+1;

  % Make Hankel matrices Y and U
  Y=zeros(twoi,j); U=zeros(twoi,j);
  for k=1:2*i
    Y(k,:)=y(k:k+j-1); U(k,:)=u(k:k+j-1);
  end
  % Compute the R factor
  R = triu(qr([U;Y]'))'; % R factor
  R = R(1:4*i,1:4*i); 	 % Truncate

  % ------------------------------------------------------------------
  % STEP 1: Calculate oblique and orthogonal projections
  % ------------------------------------------------------------------
  Rf = R(3*i+1:4*i,:);              % Future outputs
  Rp = [R(1:1*i,:);R(2*i+1:3*i,:)]; % Past inputs and outputs
  Ru  = R(1*i+1:2*i,1:twoi); 	      % Future inputs
  % Perpendicular future outputs 
  Rfp = [Rf(:,1:twoi) - (Rf(:,1:twoi)/Ru)*Ru,Rf(:,twoi+1:4*i)]; 
  % Perpendicular past inputs and outputs
  Rpp = [Rp(:,1:twoi) - (Rp(:,1:twoi)/Ru)*Ru,Rp(:,twoi+1:4*i)]; 

  % The oblique projection is computed as (6.1) in VODM, page 166.
  % obl/Ufp = Yf/Ufp * pinv(Wp/Ufp) * (Wp/Ufp)
  % The extra projection on Ufp (Uf perpendicular) tends to give 
  % better numerical conditioning (see algo on VODM page 131)

  % Funny rank check (SVD takes too long)
  % This check is needed to avoid rank deficiency warnings
  if (norm(Rpp(:,3*i-2:3*i),'fro')) < 1e-10
    Ob = (Rfp*pinv(Rpp')')*Rp; 	% Oblique projection
  else
    Ob = (Rfp/Rpp)*Rp;
  end

  % ------------------------------------------------------------------
  % STEP 2: Compute weighted oblique projection and its SVD
  %         Extra projection of Ob on Uf perpendicular
  % ------------------------------------------------------------------
  WOW = [Ob(:,1:twoi) - (Ob(:,1:twoi)/Ru)*Ru,Ob(:,twoi+1:4*i)];
  [U,S,~] = svd(WOW);
  ss = diag(S);

  % ------------------------------------------------------------------
  % STEP 3: Partitioning U into U1 and U2 (the latter is not used)
  % ------------------------------------------------------------------
  U1 = U(:,1:n); % Determine U1

  % ------------------------------------------------------------------
  % STEP 4: Determine gam = Gamma(i) and gamm = Gamma(i-1) 
  % ------------------------------------------------------------------
  gam  = U1*diag(sqrt(ss(1:n)));
  gamm = gam(1:(i-1),:);
  gam_inv  = pinv(gam); 			% Pseudo inverse of gam
  gamm_inv = pinv(gamm); 			% Pseudo inverse of gamm

  % ------------------------------------------------------------------
  % STEP 5: Determine A matrix (also C, which is not used) 
  % ------------------------------------------------------------------
  Rhs = [gam_inv*R(3*i+1:4*i,1:3*i),zeros(n,1); R(i+1:twoi,1:3*i+1)];
  Lhs = [gamm_inv*R(3*i+1+1:4*i,1:3*i+1); R(3*i+1:3*i+1,1:3*i+1)];
  sol = Lhs/Rhs;    % Solve least squares for [A;C]
  A = sol(1:n,1:n); % Extract A
return
