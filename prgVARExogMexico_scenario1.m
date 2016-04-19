%% prgVARExogMexico_scenario2.m
%
%  Author: Jorge A. Chan-Lau
%  Current version: 04/11/2016
%  Previous version: 04/07/2016
%  
%  Requires:
%       Matlab Econometrics Toolbox (different than Le Sage toolbox)
%       mvrandn.m file, Zdravko Botev
%       
%  Compatible versions: 2014b, 2015a
%  Builds on prbMexicoVAR.m, which also includes full VAR model
%
%
%  Use: Estimates VAR for domestic macro-financial variables in Mexico
%
%  Data input: HaverData.xls
%  Sheets: US variables  
%          Mexico variables
%
%  Two step VAR estimation
%
%       Unrestricted VAR for US variables
%       Unrestricted VAR for Mexico variables, US variables as exogenous
%
%  Other tricks:
%
%       plotyy
%
%% Variable description
%  --------------------
%
%  US variables
%  ------------
%
%  1.    UST 3-month bid yield, secondary market, in percent
%  2.    UST 3-month bid yield, constant maturity, in percent
%  3.    UST 10 year yield, constant maturity, in percent
%  4.    CBOE VIX index, level
%  5.    Dow Jones 30, average price, average close
%  6.    Dow Jones 30, average price, close end of period
%  7.    Industrial production index, 2012 = 100
%  8.    West Texas Intermediate, Cushing, spot price FOB avg, USD/barrel
%  9.    WTI, Cushing, spot price FOB, USD/barrel, eop
%  10.   WTI, domestic spot price, USD/barrel (CME)
%
%  Mexico variables
%  ----------------
% 
%  1.    Industrial production, 2008 = 100, in level
%  2.    CETE 28 day, rate
%  3.    MXN-USD, end-of-period, level
%  4.    MXN-USD, new peso to US$
%  5.    Consumer Price Index (CPI), level
%  6.    Unemployment rate, in percent SA, average
%  7.    Unemployment rate, in percent SA, eop
%  8.    IPC, Stock price index, level, average
%  9.    IPC, Stock price index, level, eop
%  10.   Commercial bank credit, local currency, millions, eop
%  11.   Commercial bank total performing loans, loc cur, eop, millions,
%  12.   Comm bank total performing direct loans, loc cur, eop, millions
%  13.   Comm bank tot peforming consumption loans, loc cur, eop,  millions
%  14.   Comm bank tot peforming mortgage loans, loc cur, eop,  millions
%  15.   Comm bank tot peforming firms self-empl loans, loc cur, eop,  millions
%  16.   Comm bank tot peforming non-bank fin int, loc cur, eop,  millions
%  17.   IGAE, Indicen General de Actividad Economica, level
%
%% Section I: Read Data
clc; clear all;
filename = 'HaverData.xls';
sheetname = 'US variables';

write_results = false;
plot_figures = false;

% Read US data
[A B]=xlsread(filename, sheetname);

% Variables in international VAR
%
%   2.  UST 3-month yield constant maturity
%   3.  UST 10=year yield constant maturity
%   4.  CBOE VIX index
%   6.  Dow Jones index, end period
%   7.  Industrial production, end period
%   10. WTI domestic price, CME
%
% Ordering:
%   IP, WTI, 3month, 10 year, VIX, DJ
%  

Y_us = A(:,[7 10 2 3 4 6]);
date_str = B(3:end,2);
date_str = char(date_str);
date_num = datenum(date_str);

% create data table for US data
dataUS = table(date_str,A(:,7),A(:,10),A(:,2),A(:,3),A(:,4),A(:,6), ...
    'VariableNames',{'Date','IndProd','WTI','yld3mo','yld10yr','VIX','DJ'});

% log differences
dY_us01 = log(Y_us(2:end,:)./Y_us(1:end-1,:));      % monthly
dY_us03 = log(Y_us(4:end,:)./Y_us(1:end-3,:));      % quarterly
dY_us12 = log(Y_us(13:end,:)./Y_us(1:end-12,:));    % annual

% no overlapping log differences
dY_us03novlp = dY_us03(1:3:end,:);
dY_us12novlp = dY_us12(1:12:end,:);

% Read Mexican data

sheetname = 'Mexico variables';
[A B]=xlsread(filename, sheetname);

% Variables in international VAR
%
%   3.  MXN-USD exchange rate, loc cur/USD
%   2.  Cete 28 rate, percent
%   9.  IPC Stock market index
%   23. IGAE, Indice General de Actividad Economica
%   10. Commercial loans, total, loc cur
%   5.  CPI, consumer price index, level
%   7.  Unemployment, rate

Y_mx = A(:,[3 2 9 23 10 5 7]);

% create data table for Mexico
dataMX = table(date_str,Y_mx(:,1),Y_mx(:,2),Y_mx(:,3),Y_mx(:,4),...
    Y_mx(:,5),Y_mx(:,6),Y_mx(:,7), ...
    'VariableNames',{'Date','FX','Cete28','IPC','IGAE','Credit', ...
    'CPI','UEMP'});

dY_mx01 = log(Y_mx(2:end,:)./Y_mx(1:end-1,:));      % monthly
dY_mx03 = log(Y_mx(4:end,:)./Y_mx(1:end-3,:));      % quarterly
dY_mx12 = log(Y_mx(13:end,:)./Y_mx(1:end-12,:));    % annual

% no overlapping log differences
dY_mx03novlp = dY_mx03(1:3:end,:);
dY_mx12novlp = dY_mx12(1:12:end,:);

%% Section II: US VAR model

nUS = size(dY_us01,2);
nAR = 4;
cons = true;
Y = dY_us01;
idxUSnotNaN = find(~isnan(sum(Y,2)));   % used only valid observations
Y = Y(idxUSnotNaN,:);

%  model setup
mdlUS = vgxset('n',nUS,'nAR',nAR,'Constant',cons);

[VARmdlUS, VARmdlUSStdErrors, logLmdlUS, W] = ...
    vgxvarx(mdlUS,Y(5:end,:),[],Y(1:4,:));

%auxspec = vgxvarx(mdlUS, Y(5:end,:),[],Y(1:4,:));
%T=36;
%[Forecast, ForecastCov] = vgxpred(VARmdlUS,T,[], Y(end-3:end,:));
%vgxplot(VARmdlUS, Y(end-12:end,:), Forecast, ForecastCov);


%% Section II: Mexico VAR model 2, exogenous US variables

nMX = size(dY_mx01,2);
nAR = 4;
cons = true;
Y = dY_mx01;
idxMXnotNaN = find(~isnan(sum(Y,2)));   % index valid obs
Y = Y(idxMXnotNaN,:);                   % use only valid obss

% Lagged values of X
lagvector = [0 1 2 3 4];
X = lagmatrix(dY_us01, lagvector);
X = X(idxMXnotNaN,:);

% Create design matrix, cell structure
Xexpand = kron(X,eye(nMX));
Xcell = mat2cell(Xexpand, nMX*ones(size(X,1),1));
nX = size(Xexpand,2);

% Model setup
mdlMX2 = vgxset('n',nMX,'nAR',nAR,'nX',nX,'Constant', true);   
[EstSpec, EstStdErrors, LLF, W] = vgxvarx(mdlMX2,Y,Xcell);
VARMX2 = struct('mdlEst',EstSpec,'stdErr',EstStdErrors, ...
    'logLik',LLF,'shock',W);


%% Section IV: Baseline Scenario
%
%  US growth at 2 percent annualy
%   - 
%
%  Mexico shocks consistent with Banco de Mexico FSR
%   - 1 std for all variables except
%   - unemployment: 90 percent of maximum during GFC
%   - IGAE, 2 std dev shock
%  
%  Variables suffixes
%  Xorig, Yorig:        historical time series
%  Xbase, Ybase:        months until end of 2016 (9 months)
%  Xextend, Yextend     [Series orig; Series extend]
%  Xfcst, Yfcst:        24 months after end of 2016
%  Xtot, Ytot:          [Xorig; Xbase; Xfcst]


% US 2016 baseline
% Last observed data: March 2006
% Length of 2016 baseline: 9 months
% ---------------------------------

% General parameters
base_periods = 9;   % end of 2017
T = 24;              % remainder periods until 2018

% Historical data, fill NaN data with last value
Xorig = table2array(dataUS(:,2:end));
Xorig(end,1) = Xorig(end-1,1);          % assume zero growth IP in 3/2016

% Baseline 2016
Xbase = zeros(base_periods,size(Xorig,2));

% X(:,1) Industrial production 
% growth by 2 percent per annum during 2016
IndProdGrowthAnnum = 0.02;         % 0.02 if growing
IndProdGrowthMonth = (1+IndProdGrowthAnnum)^(1/12);
IndProdGrowthMonth = repmat(IndProdGrowthMonth, base_periods,1);
IndProdGrowthMonth = cumprod(IndProdGrowthMonth);
Xbase(:,1) = Xorig(end,1)*IndProdGrowthMonth;

% X(:,2) WTI
% U.S. Energy Information Administration Forecast
% 2016 34.04
% 2017 40.09
WTI2016 = 34.04;
Xbase(end,2) = WTI2016;
step = (Xbase(end,2) - Xorig(end,2))/base_periods;
Xbase(:,2) = cumsum(repmat(step,base_periods,1)) + Xorig(end,2);

% X(:,3) US 3-month yield
% 2016 100 bps hike
shock = 1;    % 1 if 100 bps hike
Xbase(end,3) = shock + Xorig(end,3);
step = (Xbase(end,3) - Xorig(end,3))/base_periods;
Xbase(:,3) = cumsum(repmat(step,base_periods,1))+Xorig(end,3);

% X(:,4) US 10-year yield
% 2016 200 bps hike
shock = 2;   % 2 if 200 bps hike
step = shock/base_periods;
Xbase(:,4) = cumsum(repmat(step,base_periods,1))+Xorig(end,4);

% X(:,5) VIX
% 2016 Up by 1 standard deviation
stdVIX = std(dataUS.VIX);
shockfactor = 1;
shock = shockfactor*stdVIX;
step = shock/base_periods;
Xbase(:,5) = cumsum(repmat(step,base_periods,1))+Xorig(end,5);

% X(:,6) Dow Jones
% 2016 Down by 1 standard deviation, growth rate
stdDow = -std(dY_us12(:,6));
shockfactor = 1;
shock = shockfactor*stdDow;
step = (1+shock)^(1/base_periods);
Xbase(:,6) = cumprod(repmat(step,base_periods,1))*Xorig(end,6);

% Extend historical series by 2016 baseline
% calculate 1-month growth rates
Xextend = [Xorig; Xbase];
dXextend = diff(log(Xextend));

% T-month forecast for the remaining values
[fcst, fcastCov] = vgxpred(VARmdlUS,T,[], dXextend(end-3:end,:));
% vgxplot(VARmdlUS, dXextend(end-12:end,:), Forecast, ForecastCov);
Xfcst=exp(cumsum(fcst));
Xfcst=bsxfun(@times,Xextend(end,:),Xfcst);
Xtot = [Xextend; Xfcst];
dXtot = diff(log(Xtot));

% Mexico 2016 baseline
% Last observed data: March 2006
% Length of 2016 baseline: 9 months
% ---------------------------------

% Complete NaN
% Variable 4:   IGAE (last two periods incomplete)
% Variable 5:   credit
% Variable 6:   CPI
% Variable 7:   Unemployment

Yorig = table2array(dataMX(:,2:end));

Yorig(end-1,4) = Yorig(end-2,4);
Yorig(end,4) = Yorig(end-1,4);
Yorig(end,5) = Yorig(end-1,5);
Yorig(end,6) = Yorig(end-1,6);
Yorig(end,7) = Yorig(end-1,7);

% Mexico 2016 baseline
Ybase=zeros(base_periods, size(Yorig,2));
stdshock = 1.0;  % standard deviation shock, all except IGAE
stdIGAE  = 2.5;  % standard deviation shock, IGAE

mx_var_array = [1 2 3 5 6];
dY_array = dY_mx12;

% negative shocks to IPC, Credit, and CPI
negative_shock = [3 5];

for i=1:length(mx_var_array),
    idx = mx_var_array(i);
    shock = nanstd(dY_array(:,idx))*stdshock;   
    if ismember(idx,negative_shock),
        shock = -shock;
    end
    step = (1+shock)^(1/base_periods);
    step = cumprod(repmat(step,base_periods,1));
    Ybase(:,idx) = Yorig(end,idx)*step;
end,

% Unemployment
step = (0.9*nanmax(Yorig(:,7))-Yorig(end,7))/base_periods;
Ybase(:,7) = cumsum(repmat(step,base_periods,1))+Yorig(end,7);

% IGAE
shock = -nanstd(dY_array(:,4))*stdIGAE;   
step = (1+shock)^(1/base_periods);
step = cumprod(repmat(step,base_periods,1));
Ybase(:,4) = Yorig(end,4)*step;

Yextend = [Yorig; Ybase];           % Add base_periods
dYextend = diff(log(Yextend));

% 24-month forecast for remainder of forecast horizon
% ----------------------------------------------------

% Create design matrix, exogenous regressors, cell structure

idxMXvalid = find(~isnan(sum(dYextend,2)));          % index valid obs
dY = dYextend(idxMXvalid,:);                         % only Y valid obs
idxMXvalid = [idxMXvalid; (1:1:T)'+idxMXvalid(end)]; % extend to Xtot

lagvector = [0 1 2 3 4];                        % lags used in regression
X = lagmatrix(dXtot, lagvector);             % lagged matrix
X = X(idxMXvalid,:);                            % only valid obs
Xexpand = kron(X,eye(nMX));                     % cell structure
Xcell = mat2cell(Xexpand, nMX*ones(size(X,1),1));
nX = size(Xexpand,2);

[Yfcst, YfcstCov] = vgxpred(VARMX2.mdlEst,T,Xcell, dYextend(end-3:end,:));

Yfcst=exp(cumsum(Yfcst));
Yfcst=bsxfun(@times,Yextend(end,:),Yfcst);
Ytot = [Yextend; Yfcst];
dYtot = diff(log(Ytot));
dYtot12 = log(Ytot(13:end,:)./Ytot(1:end-12,:));

% Estimate innovations to domestic variables
% in the baseline 2016 scenario

idxMXtot = find(~isnan(sum(dYtot,2)));
[W, logL] = vgxinfer(VARMX2.mdlEst,dYtot(idxMXtot,:), Xcell);
Wbase = W(end-T-base_periods+1:end-T,:);
Whist = W(1:end-T-base_periods,:);
Wstd  = nanstd(Whist);

% Generate innovations
rng(10);     % for reproducibility
nsim = 2000;

% Base period, remainder of 2016, using truncated normal random generator
shift_factor = 0.20;
%Wbase_shocks_std = bsxfun(@ldivide,Wstd,Wbase);  % standard normal shocks
Wbase_shocks_minus = (1-shift_factor)*Wbase;    
Wbase_shocks_plus  = (1+shift_factor)*Wbase;
Wbase_lowbnd = min(Wbase_shocks_minus, Wbase_shocks_plus);
Wbase_uppbnd = max(Wbase_shocks_minus, Wbase_shocks_plus);
%Wbase_lowbnd = bsxfun(@times,Wstd, Wbase_low_bound);
%Wbase_uppbnd = bsxfun(@times,Wstd, Wbase_upp_bound);

Q = VARMX2.mdlEst.Q; 
base_sim = zeros(base_periods,nsim*nMX);
for i=1:base_periods,
    % drawing from truncated random distribution
    aux=mvrandn(Wbase_lowbnd(i,:)',Wbase_uppbnd(i,:)',Q,nsim)';
    aux=reshape(aux',1,numel(aux));
    base_sim(i,:) = aux;
end,

% Forecast period, after 2016
lb = zeros(nMX,1) - Inf;
ub = zeros(nMX,1) + Inf;
fcst_sim = zeros(T,nsim*nMX);
for i=1:T,
    % drawing from truncated random distribution
    aux=mvrandn(lb,ub,Q,nsim)';
    aux=reshape(aux',1,numel(aux));
    fcst_sim(i,:) = aux;
end,

% innovations
% (base_periods + T) x (nMX x nsim)
Wsim = [base_sim; fcst_sim];
Wsim3D = reshape(Wsim,base_periods+T,nMX,nsim);

[Ysim, logL] = vgxproc(VARMX2.mdlEst,Wsim3D,...
    Xcell(end-base_periods-T+1:end,:), ...
    dYtot(end-base_periods-T-4+1:end,:), ...
    Whist(end-4+1:end,:));

%% Section IV - create scenario dates

months = 1:1:12;
years = 2016:1:2018;
nmonths = 12;
nyears = 3;

% Create date series

day = zeros(nmonths*nyears,1);
scenario_dates = cell(nmonths*nyears,1);
fmtDate = 'mm/dd/yyyy';

for i=1:nyears,
    for j=1:nmonths,
        thisCounter = (i-1)*12 + j;
        day(thisCounter) = eomday(years(i),months(j));
        scenario_dates{thisCounter}= ...
            datestr([years(i),months(j),day(thisCounter),0,0,0], fmtDate);
    end, 
end,
scenario_dates = char(scenario_dates);
date_sim = [date_str; scenario_dates(4:end,:)];

%% Simulation paths, write results
%
if write_results, 
    
    % Adverse scenarios, in levels
        
    outfile = 'MXAdverse.xlsx';
    sheets = cell(nMX,1);
    sheets{1} = 'Exchange Rate';
    sheets{2} = 'Cete 28';
    sheets{3} = 'IPC';
    sheets{4} = 'IGAE';
    sheets{5} = 'Credit';
    sheets{6} = 'CPI';
    sheets{7} = 'Unemployment';

    % Name of scenarios: Sim 1, Sim 2, etc
    strname = {'Sim '};
    simul_name = cellstr(strcat(strname,'1'));
    for i = 2:nsim,
        next_simul = cellstr(strcat(strname,num2str(i)));
        simul_name = [simul_name, next_simul];
    end, 

    n_dates = size(Ysim,1);    
    horizon = 36;
    
    b= {'Date'};
    for i = 1:nMX,
        idx_var = i;
        var_mx = squeeze(Ysim(:,idx_var,:));
        h = exp(cumsum(var_mx))*Yorig(end,idx_var);
        h = [repmat(Yorig(:,idx_var), 1, size(h,2)); h];
        xlswrite(outfile,h(end-horizon+1:end,:),sheets{i},'B2');  
        xlswrite(outfile,cellstr(date_sim(end-horizon+1:end,:)),...
            sheets{i},'A2');    
        xlswrite(outfile,simul_name, sheets{i},'B1');
        xlswrite(outfile,b, sheets{i},'A1');
    end,

    % Adverse scenarios, growth rates

    outfile = 'MXAdverseGrowth.xlsx';
    for i = 1:nMX,
        idx_var = i;
        var_mx = squeeze(Ysim(:,idx_var,:));
        h = exp(cumsum(var_mx))*Yorig(end,idx_var);
        % variable in levels
        h = [repmat(Yorig(:,idx_var), 1, size(h,2)); h];
        % 12-month growth rate
        hdif = log(h(13:end,:)./h(1:end-12,:));   
        xlswrite(outfile,hdif(end-horizon+1:end,:),sheets{i},'B2');  
        xlswrite(outfile,cellstr(date_sim(end-horizon+1:end,:)), ...
            sheets{i},'A2');    
        xlswrite(outfile,simul_name, sheets{i},'B1');
        xlswrite(outfile,b, sheets{i},'A1');
    end,

end, % end write_file

%% Table with min, max, and average shocks

proj_level_max =[];
proj_level_min =[];
proj_level_mean = [];
proj_growth_max =[];
proj_growth_min =[];
proj_growth_mean = [];

for i = 1:nMX,
        idx_var = i;
        var_mx = squeeze(Ysim(:,idx_var,:));
        h = exp(cumsum(var_mx))*Yorig(end,idx_var);
        % variable in levels
        h = [repmat(Yorig(:,idx_var), 1, size(h,2)); h];
        hmax = max(h,[],2);
        hmin = min(h,[],2);
        hmean = mean(h,2);
        
        proj_level_max = [proj_level_max, hmax];
        proj_level_min = [proj_level_max, hmin];
        proj_level_mean = [proj_level_max, hmean];
        
        % 12-month growth rate
        hdif = log(h(13:end,:)./h(1:end-12,:));   
        hdifmax = max(hdif,[],2);
        hdifmin = min(hdif,[],2);
        hdifmean = mean(hdif,2);
        
        proj_growth_max = [proj_growth_max, hdifmax];
        proj_growth_min = [proj_growth_max, hdifmin];
        proj_growth_mean = [proj_growth_max, hdifmean];
        
        
end, 

horizon = 36;
end_period = size(proj_level_max,1);
idx2016 = (end_period-horizon+1):1:(end_period-horizon+12);
idx2017 = (end_period-horizon+12+1):1:(end_period-horizon+24);
idx2018 = (end_period-horizon+24+1):1:(end_period);

table_results = zeros(nMX,7);

level_set_high=[1,2,7];
level_set_low =[3];
growth_set_low = [4,5];
growth_set_high= [6];

for i = 1:nMX,
    idx_var = i;
    if ismember(idx_var,level_set_high),
        table_results(i,1) = proj_level_max(end_period-horizon,idx_var);
        table_results(i,2) = max(proj_level_max(idx2016,idx_var));
        table_results(i,3) = mean(proj_level_mean(idx2016,idx_var));
        table_results(i,4) = max(proj_level_max(idx2017,idx_var));
        table_results(i,5) = mean(proj_level_mean(idx2017,idx_var));
        table_results(i,6) = max(proj_level_max(idx2018,idx_var));
        table_results(i,7) = mean(proj_level_mean(idx2018,idx_var));       
    end,
    if ismember(idx_var,level_set_low),
        table_results(i,1) = proj_level_max(end_period-horizon,idx_var);
        table_results(i,2) = min(proj_level_min(idx2016,idx_var));
        table_results(i,3) = mean(proj_level_mean(idx2016,idx_var));
        table_results(i,4) = min(proj_level_min(idx2017,idx_var));
        table_results(i,5) = mean(proj_level_mean(idx2017,idx_var));
        table_results(i,6) = min(proj_level_min(idx2018,idx_var));
        table_results(i,7) = mean(proj_level_mean(idx2018,idx_var));       
    end
    if ismember(idx_var,growth_set_high),
        table_results(i,1) = proj_growth_max(end_period-horizon-12,idx_var);
        table_results(i,2) = max(proj_growth_max(idx2016-12,idx_var));
        table_results(i,3) = mean(proj_growth_mean(idx2016-12,idx_var));
        table_results(i,4) = max(proj_growth_max(idx2017-12,idx_var));
        table_results(i,5) = mean(proj_growth_mean(idx2017-12,idx_var));
        table_results(i,6) = max(proj_growth_max(idx2018-12,idx_var));
        table_results(i,7) = mean(proj_growth_mean(idx2018-12,idx_var));       
    end,
    if ismember(idx_var,growth_set_low),
        table_results(i,1) = proj_growth_max(end_period-horizon-12,idx_var);
        table_results(i,2) = min(proj_growth_min(idx2016-12,idx_var));
        table_results(i,3) = mean(proj_growth_mean(idx2016-12,idx_var));
        table_results(i,4) = min(proj_growth_min(idx2017-12,idx_var));
        table_results(i,5) = mean(proj_growth_mean(idx2017-12,idx_var));
        table_results(i,6) = min(proj_growth_min(idx2018-12,idx_var));
        table_results(i,7) = mean(proj_growth_mean(idx2018-12,idx_var));       
    end
 
end


%% Plot figures

if plot_figures, 
    
    date_sim = datenum(date_sim); % create dates
    date_12m = date_sim(13:end,:);
    end_orig =size(Yorig,1)-12;

    % Label for figures, growth rates
    ylab = cell(7,1);
    ylab{1} = 'MXN-USD exchange rate';
    ylab{2} = 'Cete 28-day rate, in percent';
    ylab{3} = 'IPC, annual growth rate, in percent';
    ylab{4} = 'IGAE, annual growth rate, in percent';
    ylab{5} = 'Credit, annual growth rate, in percent';
    ylab{6} = 'CPI, annual growth rate, in percent';
    ylab{7} = 'Unemployment rate,in percent';
    scen = {'Adverse scenario: '};

    % Growth rates

    for i = 1:nMX,
        idx_var = i;
        var_mx = squeeze(Ysim(:,idx_var,:));
        h = exp(cumsum(var_mx))*Yorig(end,idx_var);
        % variable in levels
        h = [repmat(Yorig(:,idx_var), 1, size(h,2)); h];
        % 12-month growth rate
        hdif = log(h(13:end,:)./h(1:end-12,:));   
        hdifmax = max(hdif,[],2);
        hdifmin = min(hdif,[],2);
        hdifmean = mean(hdif,2);

        figure;
        plot(date_12m(1:end_orig),hdif(1:end_orig),'-k'); hold on;
        plot(date_12m(end_orig+1:end),...
            [hdifmax(end_orig+1:end), hdifmin(end_orig+1:end)], ...
            'r','LineWidth',1.5); 
        hold on;
        plot(date_12m(end_orig+1:end),hdifmean(end_orig+1:end), '-b',...
            'LineWidth',2); hold off;    
        set(gca,'Xtick',[date_12m([1:12:length(date_12m)])]); 
        xlim([date_12m(1); date_12m(end)]); 
        ylim([1.1*min(hdifmin) 1.1*max(hdifmax)]);
        rline = refline([0 0]);
        line([date_12m(end_orig) date_12m(end_orig)],...
            [1.1*min(hdifmin) 1.1*max(hdifmax)], ...
            'Color','k', 'LineStyle','--');    
        datetick('x', 'mmmYY','keepticks','keeplimits');             
        ylabel(ylab{i}, 'FontSize',14);
        title(strcat(scen, ylab{i}), 'FontSize', 16, 'FontWeight','normal');
        ax = gca;
        ax.PlotBoxAspectRatio = [4 1 1];

    end,

    % Label for figures, levels
    ylab{1} = 'MXN-USD exchange rate';
    ylab{2} = 'Cete 28-day rate, in percent';
    ylab{3} = 'IPC,in levels';
    ylab{4} = 'IGAE, in levels';
    ylab{5} = 'Credit, in millions of local currency';
    ylab{6} = 'CPI';
    ylab{7} = 'Unemployment rate,in percent';
    scen = {'Adverse scenario: '};

    % Levels
    date_12m = date_sim;
    end_orig =size(Yorig,1);
    
    for i = 1:nMX,
        idx_var = i;
        var_mx = squeeze(Ysim(:,idx_var,:));
        h = exp(cumsum(var_mx))*Yorig(end,idx_var);
        % variable in levels
        h = [repmat(Yorig(:,idx_var), 1, size(h,2)); h];

        hmax = max(h,[],2);
        hmin = min(h,[],2);
        hmean = mean(h,2);

        figure;
        plot(date_12m(1:end_orig),h(1:end_orig),'-k'); hold on;
        plot(date_12m(end_orig+1:end),...
            [hmax(end_orig+1:end), hmin(end_orig+1:end)], ...
            'r','LineWidth',1.5); 
        hold on;
        plot(date_12m(end_orig+1:end),hmean(end_orig+1:end), '-b',...
            'LineWidth',2); hold off;    
        set(gca,'Xtick',[date_12m([1:12:length(date_12m)])]); 
        xlim([date_12m(1); date_12m(end)]); 
        ylim([0.9*min(hmin) 1.1*max(hmax)]);
        line([date_12m(end_orig) date_12m(end_orig)],...
            [0.9*min(hmin) 1.1*max(hmax)], ...
            'Color','k', 'LineStyle','--');    
        datetick('x', 'mmmYY','keepticks','keeplimits');         
        ylabel(ylab{i}, 'FontSize',14);
        title(strcat(scen, ylab{i}), 'FontSize', 16, 'FontWeight','normal');
        ax = gca;
        ax.PlotBoxAspectRatio = [4 1 1];


    end,

end, % end if plot_figures

%% Average paths

Yavg = [];
Yavggrowth =[];
for i = 1:nMX,

    idx_var = i;
    var_mx = squeeze(Ysim(:,idx_var,:));
    h = exp(cumsum(var_mx))*Yorig(end,idx_var);
    % variable in levels
    h = [repmat(Yorig(:,idx_var), 1, size(h,2)); h];

    hmax = max(h,[],2);
    hmin = min(h,[],2);
    hmean = mean(h,2);
  
    Yavg=[Yavg, hmean];

end, 

Yavggrowth = [];
for i = 1:nMX,
    idx_var = i;
    var_mx = squeeze(Ysim(:,idx_var,:));
    h = exp(cumsum(var_mx))*Yorig(end,idx_var);
    % variable in levels
    h = [repmat(Yorig(:,idx_var), 1, size(h,2)); h];
    % 12-month growth rate
    hdif = log(h(13:end,:)./h(1:end-12,:));   
    hdifmax = max(hdif,[],2);
    hdifmin = min(hdif,[],2);
    hdifmean = mean(hdif,2);
    Yavggrowth = [Yavggrowth, hdifmean];
end, 

%% Loose ends
quartile_vector = [0.05 0.25 0.50 0.75 0.95];
h = squeeze(Ysim(:,4,:));
h = exp(cumsum(h));
hq = quantile(h, [0.05, 0.25, 0.50, 0.75, 0.95],2);

% keep only simulations below a certain percentile
idx_pctl = 2;   % 25 percentile
numbelow = 20; % accept if at least these number of observations are below

h_shrt = h(base_periods + 1:T+base_periods,:);
hq_shrt = h(base_periods + 1:T+base_periods,:);

idx = bsxfun(@gt, hq_shrt(:,idx_pctl), h_shrt);
idx = find(sum(idx,1)>=T);

% Plot individual series
idx_series = 4;  %IGAE
h = exp(cumsum( squeeze(Ysim(:,idx_series,:))))*Yorig(end,idx_series);
h = [repmat(Yorig(:,idx_series), 1, size(h,2)); h];
xaxis = 1:1:(size(Yorig,1)+base_periods+T);

hmax = max(h,[],2);     % upper bound of projections
hmin = min(h,[],2);     % lower bound of projections
hmean = mean(h,2);   % average path
hmaxgrowth = log(hmax(13:end,:)./hmax(1:end-12,:));
hmingrowth = log(hmin(13:end,:)./hmin(1:end-12,:));
hmeangrowth = log(hmean(13:end,:)./hmean(1:end-12,:));
hq=quantile(h,quartile_vector,2);
plot([hmin hmax hmean]);
plot(hq);
plot([hmingrowth hmaxgrowth hmeangrowth]);
