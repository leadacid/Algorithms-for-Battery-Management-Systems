% --------------------------------------------------------------------
% function processOCV
%
% Technical note: PROCESSOCV assumes that specific Arbin test scripts
% have been executed to generate the input files. "makeMATfiles.m" 
% converts the raw Excel data files into "MAT" format, where the MAT 
% files have fields for time, step, current, voltage, chgAh, and disAh
% for each script run.
%
% The results from four scripts are required at every temperature.  
% The steps in each script file are assumed to be:
%   Script 1 (thermal chamber set to test temperature):
%     Step 1: Rest @ 100% SOC to acclimatize to test temperature
%     Step 2: Discharge @ low rate (ca. C/30) to min voltage
%     Step 3: Rest ca. 0%
%   Script 2 (thermal chamber set to 25 degC):
%     Step 1: Rest ca. 0% SOC to acclimatize to 25 degC
%     Step 2: Discharge to min voltage (ca. C/3)
%     Step 3: Rest
%     Step 4: Constant voltage at vmin until current small (ca. C/30)
%     Steps 5-7: Dither around vmin
%     Step 8: Rest
%     Step 9: Constant voltage at vmin for 15 min
%     Step 10: Rest
%   Script 3 (thermal chamber set to test temperature):
%     Step 1: Rest at 0% SOC to acclimatize to test temperature
%     Step 2: Charge @ low rate (ca. C/30) to max voltage
%     Step 3: Rest
%   Script 4 (thermal chamber set to 25 degC):
%     Step 1: Rest ca. 100% SOC to acclimatize to 25 degC
%     Step 2: Charge to max voltage (ca. C/3)
%     Step 3: Rest
%     Step 4: Constant voltage at vmax until current small (ca. C/30)
%     Steps 5-7: Dither around vmax
%     Step 8: Rest
%     Step 9: Constant voltage at vmax for 15 min
%     Step 10: Rest
% 
% All other steps (if present) are ignored by PROCESSOCV.  The time 
% step between data samples is not critical since the Arbin integrates
% ampere-hours to produce the two Ah columns, and this is what is 
% necessary to generate the OCV curves.  The rest steps must 
% contain at least one data point each.

function model=processOCV(data,cellID,minV,maxV,savePlots)
  filetemps = [data.temp]; filetemps = filetemps(:);
  numtemps = length(filetemps); 
  
  ind25 = find(filetemps == 25);
  if isempty(ind25),
    error('Must have a test at 25degC');
  end
  not25 = find(filetemps ~= 25);

  % ------------------------------------------------------------------
  % Process 25 degC data to find raw OCV relationship and eta25
  % ------------------------------------------------------------------
  SOC = 0:0.005:1; % output SOC points for this step
  filedata = zeros([0 length(data)]);
  eta = zeros(size(filetemps)); % coulombic efficiency
  Q   = zeros(size(filetemps)); % apparent total capacity
  data25 = data(ind25);

  totDisAh = data25.script1.disAh(end) + ...
             data25.script2.disAh(end) + ...
             data25.script3.disAh(end) + ...
             data25.script4.disAh(end);
  totChgAh = data25.script1.chgAh(end) + ...
             data25.script2.chgAh(end) + ...
             data25.script3.chgAh(end) + ...
             data25.script4.chgAh(end);
  eta25 = totDisAh/totChgAh; eta(ind25) = eta25;
  data25.script1.chgAh = data25.script1.chgAh*eta25;
  data25.script2.chgAh = data25.script2.chgAh*eta25;
  data25.script3.chgAh = data25.script3.chgAh*eta25;
  data25.script4.chgAh = data25.script4.chgAh*eta25;

  Q25 = data25.script1.disAh(end) + data25.script2.disAh(end) - ...
        data25.script1.chgAh(end) - data25.script2.chgAh(end);
  Q(ind25) = Q25;
  indD  = find(data25.script1.step == 2); % slow discharge
  IR1Da = data25.script1.voltage(indD(1)-1) - ...
          data25.script1.voltage(indD(1));
  IR2Da = data25.script1.voltage(indD(end)+1) - ...
          data25.script1.voltage(indD(end));
  indC  = find(data25.script3.step == 2);
  IR1Ca = data25.script3.voltage(indC(1)) - ...
          data25.script3.voltage(indC(1)-1);
  IR2Ca = data25.script3.voltage(indC(end)) - ...
          data25.script3.voltage(indC(end)+1);
  IR1D = min(IR1Da,2*IR2Ca); IR2D = min(IR2Da,2*IR1Ca);
  IR1C = min(IR1Ca,2*IR2Da); IR2C = min(IR2Ca,2*IR1Da);
  
  blend = (0:length(indD)-1)/(length(indD)-1);
  IRblend = IR1D + (IR2D-IR1D)*blend(:);
  disV = data25.script1.voltage(indD) + IRblend;
  disZ = 1 - data25.script1.disAh(indD)/Q25;
  disZ = disZ + (1 - disZ(1));
  filedata(ind25).disZ = disZ; 
  filedata(ind25).disV = data25.script1.voltage(indD);
  filedata(ind25).disVmod = disV;
  
  blend = (0:length(indC)-1)/(length(indC)-1);
  IRblend = IR1C + (IR2C-IR1C)*blend(:);
  chgV = data25.script3.voltage(indC) - IRblend;
  chgZ = data25.script3.chgAh(indC)/Q25;
  chgZ = chgZ - chgZ(1);
  filedata(ind25).chgZ = chgZ; 
  filedata(ind25).chgV = data25.script3.voltage(indC);
  filedata(ind25).chgVmod = chgV;

  deltaV50 = interp1(chgZ,chgV,0.5) - interp1(disZ,disV,0.5);
  ind = find(chgZ < 0.5);
  vChg = chgV(ind) - chgZ(ind)*deltaV50;
  zChg = chgZ(ind);
  ind = find(disZ > 0.5);
  vDis = flipud(disV(ind) + (1 - disZ(ind))*deltaV50);
  zDis = flipud(disZ(ind));
  filedata(ind25).rawocv = interp1([zChg; zDis],[vChg; vDis],SOC,...
                               'linear','extrap');
  
  filedata(ind25).temp = data25.temp;
  
  % ------------------------------------------------------------------
  % Process other temperatures to find raw OCV relationship and eta
  % ------------------------------------------------------------------
  for k = not25',
    data(k).script2.chgAh = data(k).script2.chgAh*eta25;
    data(k).script4.chgAh = data(k).script4.chgAh*eta25;    
    eta(k) = (data(k).script1.disAh(end) + ...
              data(k).script2.disAh(end) + ...
              data(k).script3.disAh(end) + ...
              data(k).script4.disAh(end) - ...
              data(k).script2.chgAh(end) - ...
              data(k).script4.chgAh(end))/ ...
             (data(k).script1.chgAh(end) + ...
              data(k).script3.chgAh(end));
    data(k).script1.chgAh = eta(k)*data(k).script1.chgAh;         
    data(k).script3.chgAh = eta(k)*data(k).script3.chgAh;         

    Q(k) = data(k).script1.disAh(end) + data(k).script2.disAh(end) ...
           - data(k).script1.chgAh(end) - data(k).script2.chgAh(end);
    indD = find(data(k).script1.step == 2); % slow discharge
    IR1Da = data(k).script1.voltage(indD(1)-1) - ...
            data(k).script1.voltage(indD(1));
    IR2Da = data(k).script1.voltage(indD(end)+1) - ...
            data(k).script1.voltage(indD(end));
    indC = find(data(k).script3.step == 2);
    IR1Ca = data(k).script3.voltage(indC(1)) - ...
            data(k).script3.voltage(indC(1)-1);
    IR2Ca = data(k).script3.voltage(indC(end)) - ...
            data(k).script3.voltage(indC(end)+1);
    IR1D = min(IR1Da,2*IR2Ca); IR2D = min(IR2Da,2*IR1Ca);
    IR1C = min(IR1Ca,2*IR2Da); IR2C = min(IR2Ca,2*IR1Da);

    blend = (0:length(indD)-1)/(length(indD)-1);
    IRblend = IR1D + (IR2D-IR1D)*blend(:);
    disV = data(k).script1.voltage(indD) + IRblend;
    disZ = 1 - data(k).script1.disAh(indD)/Q25;
    disZ = disZ + (1 - disZ(1));
    filedata(k).disZ = disZ; 
    filedata(k).disV = data(k).script1.voltage(indD);
    filedata(k).disVmod = disV;
    
    blend = (0:length(indC)-1)/(length(indC)-1);
    IRblend = IR1C + (IR2C-IR1C)*blend(:);
    chgV = data(k).script3.voltage(indC) - IRblend;
    chgZ = data(k).script3.chgAh(indC)/Q25;
    chgZ = chgZ - chgZ(1);
    filedata(k).chgZ = chgZ; 
    filedata(k).chgV = data(k).script3.voltage(indC);
    filedata(k).chgVmod = chgV;

    deltaV50 = interp1(chgZ,chgV,0.5) - interp1(disZ,disV,0.5);
    ind = find(chgZ < 0.5);
    vChg = chgV(ind) - chgZ(ind)*deltaV50;
    zChg = chgZ(ind);
    ind = find(disZ > 0.5);
    vDis = flipud(disV(ind) + (1 - disZ(ind))*deltaV50);
    zDis = flipud(disZ(ind));
    filedata(k).rawocv = interp1([zChg; zDis],[vChg; vDis],SOC,...
                                 'linear','extrap');
  
    filedata(k).temp = data(k).temp;
  end

  % ------------------------------------------------------------------
  % Use the SOC versus OCV data now available at each individual
  % temperature to compute an OCV0 and OCVrel relationship
  % ------------------------------------------------------------------
  % First, compile the voltages and temperatures into a single array 
  % rather than a structure
  Vraw = []; temps = []; 
  for k = 1:numtemps,
    if filedata(k).temp > 0,
      Vraw = [Vraw; filedata(k).rawocv]; %#ok<AGROW>
      temps = [temps; filedata(k).temp]; %#ok<AGROW>
    end
  end
  X = [ones(size(temps)), temps] \ Vraw;  
  model.OCV0 = X(1,:);
  model.OCVrel = X(2,:);
  model.SOC = SOC;

  % ------------------------------------------------------------------
  % Make SOC0 and SOCrel
  % ------------------------------------------------------------------
  z = -0.1:0.01:1.1; % test soc vector
  v = minV-0.01:0.01:maxV+0.01;
  socs = [];
  for T = filetemps',
    v1 = OCVfromSOCtemp(z,T,model);
    socs = [socs; interp1(v1,z,v)]; %#ok<AGROW>
  end

  SOC0 = zeros(size(v)); SOCrel = SOC0; 
  H = [ones([numtemps,1]), filetemps]; 
  for k = 1:length(v),
    X = H\socs(:,k); % fit SOC(v,T) = 1*SOC0(v) + T*SOCrel(v)
    SOC0(k) = X(1); 
    SOCrel(k) = X(2);
  end
  model.OCV = v;
  model.SOC0 = SOC0;
  model.SOCrel = SOCrel;
  
  % ------------------------------------------------------------------
  % Save data in output directory
  % ------------------------------------------------------------------
  model.OCVeta = eta;
  model.OCVQ = Q;
  model.name = cellID;
  
%   % ------------------------------------------------------------------
%   % Plot some data...
%   % ------------------------------------------------------------------
%   for k = 1:numtemps, 
%     figure(k); clf; 
%     plot(100*SOC,OCVfromSOCtemp(SOC,filedata(k).temp,model),...
%          100*SOC,filedata(k).rawocv); hold on
%     xlabel('SOC (%)'); ylabel('OCV (V)'); ylim([minV-0.1 maxV+0.1]);
%     title(sprintf('%s OCV relationship at temp = %d',...
%       cellID,filedata(k).temp)); xlim([0 100]);
%     err = filedata(k).rawocv - ...
%           OCVfromSOCtemp(SOC,filedata(k).temp,model);
%     rmserr = sqrt(mean(err.^2));
%     text(2,maxV-0.15,sprintf('RMS error = %4.1f (mV)',...
%       rmserr*1000),'fontsize',14);
%     figFormat
%     plot(100*filedata(k).disZ,filedata(k).disV,'k--','linewidth',1);
%     plot(100*filedata(k).chgZ,filedata(k).chgV,'k--','linewidth',1);
%     legend('Model prediction','Approximate OCV from data',...
%            'Raw measured data','location','southeast');
%     
%     if savePlots,
%       if ~exist('OCV_FIGURES','dir'), mkdir('OCV_FIGURES'); end
%         if filetemps(k) < 0,
%           filename = sprintf('OCV_FIGURES/%s_N%02d.eps',...
%             cellID,abs(filetemps(k)));
%         else
%           filename = sprintf('OCV_FIGURES/%s_P%02d.eps',...
%             cellID,filetemps(k));
%         end
%         print(filename,'-deps2c')
%     end
%   end

  % ------------------------------------------------------------------
  % Plot some data...
  % ------------------------------------------------------------------
  if strcmp(cellID,'E1'),
    k = find(filetemps == -15);
    figure(100); clf; 
    plot(100*filedata(k).disZ,filedata(k).disV,'-',...
         100*filedata(k).chgZ,filedata(k).chgV,'-','linewidth',1);
    hold on   
    plot(100*SOC,filedata(k).rawocv,'k--'); hold on

    xlabel('State of charge (%)'); ylabel('Voltage (V)'); ylim([minV maxV]);
    title('Voltages versus SOC');
    set(gca,'xtick',0:20:100);
    set(gca,'ytick',3:0.2:4.2);
    grid on
    bookFormatMargin;
    yt=get(gca,'YTick');
    set(gca, 'YTickLabel', num2str(yt(:), '%.1f'))  
    legend('Discharge voltage','Charge voltage','Approximate OCV',...
            'location','southeast');
    print -deps2c ../FIGURES/profileVoltage.eps

    figure(102); clf; 
    plot(100*filedata(k).disZ,filedata(k).disV,'-',...
         100*filedata(k).chgZ,filedata(k).chgV,'-','linewidth',1);
%     hold on   
%     plot(100*SOC,filedata(k).rawocv,'k--'); hold on

    xlabel('State of charge (%)'); ylabel('Voltage (V)'); ylim([minV maxV]);
    title('Voltages versus SOC');
    set(gca,'xtick',0:20:100);
    set(gca,'ytick',3:0.2:4.2);
    grid on
    notesFormat([0.1 0 0.05 0])
yt=get(gca,'YTick');
set(gca, 'YTickLabel', num2str(yt(:), '%.1f'))  
    legend('Discharge voltage','Charge voltage','location','southeast');
    print -deps2c ../FIGURES/profileVoltageTrim.eps
      
    figure(101); clf; 
    plot(100*filedata(k).disZ,filedata(k).disV,'-',...
         100*filedata(k).chgZ,filedata(k).chgV,'-','linewidth',1);
    hold on   
    plot(100*SOC,filedata(k).rawocv,'k--'); hold on

    xlabel('State of charge (%)'); ylabel('Voltage (V)'); ylim([minV maxV]);
    title('Voltages versus SOC');
    set(gca,'xtick',0:20:100);
    set(gca,'ytick',3:0.2:4.2); grid on
    notesFormat([0.1 0 0.05 0]);
yt=get(gca,'YTick');
set(gca, 'YTickLabel', num2str(yt(:), '%.1f'))  
    legend('Discharge voltage','Charge voltage','Approximate OCV',...
            'location','southeast');
    print -deps2c ~/DOSSIER/SOURCE/TEACHING/ece4710/NOTES/CH02-EquivalentCircuit/FIGURES/profileVoltage.eps
    
    figure(103); clf; 
    plot(100*filedata(k).disZ,filedata(k).disV,'-',...
         100*filedata(k).chgZ,filedata(k).chgV,'-'); hold on
    ax = gca; ax.ColorOrderIndex = 1;
    plot(100*filedata(k).disZ,filedata(k).disVmod,'--',...
         100*filedata(k).chgZ,filedata(k).chgVmod,'--');

    xlabel('State of charge (%)'); ylabel('Voltage (V)'); ylim([minV maxV]);
    title('Voltages versus SOC');
    set(gca,'xtick',0:20:100);
    set(gca,'ytick',3:0.2:4.2); grid on
    notesFormat([0.1 0 0.05 0]);
yt=get(gca,'YTick');
set(gca, 'YTickLabel', num2str(yt(:), '%.1f'))  
    legend('Discharge voltage','Charge voltage','Discharge voltage, \itR\rm removed','Charge voltage, \itR\rm removed',...
            'location','southeast');
    print -deps2c ../FIGURES/profileRremoved.eps

    figure(104); clf; 
    plot(100*filedata(k).disZ,filedata(k).disV,'-',...
         100*filedata(k).chgZ,filedata(k).chgV,'-'); hold on
    ax = gca; ax.ColorOrderIndex = 1;
    plot(100*filedata(k).disZ,filedata(k).disVmod,'--',...
         100*filedata(k).chgZ,filedata(k).chgVmod,'--');
    plot(100*SOC,filedata(k).rawocv,'k--'); hold on       

    xlabel('State of charge (%)'); ylabel('Voltage (V)'); ylim([minV maxV]);
    title('Voltages versus SOC');
    set(gca,'xtick',0:20:100);
    set(gca,'ytick',3:0.2:4.2); grid on
    notesFormat([0.1 0 0.05 0]);
yt=get(gca,'YTick');
set(gca, 'YTickLabel', num2str(yt(:), '%.1f'))  
    legend('Discharge voltage','Charge voltage','Discharge voltage, \itR\rm removed',...
      'Charge voltage, \itR\rm removed','Approximate OCV',...
            'location','southeast');
    print -deps2c ../FIGURES/approxOCV.eps
  end

end