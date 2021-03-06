%% prgVARExogMexico.m
%
%  Author: Jorge A. Chan-Lau
%  Current version: 04/08/2016
%  Previous version: 04/07/2016
%  
%  Requires:
%       Matlab Econometrics Toolbox (different than Le Sage toolbox)
%       mvrandn.m file, Zdravko Botev
%       
%  Compatible versions: 2014b, 2015a
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

auxspec = vgxvarx(mdlUS, Y(5:end,:),[],Y(1:4,:));
T=36;
[Forecast, ForecastCov] = vgxpred(VARmdlUS,T,[], Y(end-3:end,:));
vgxplot(VARmdlUS, Y(end-12:end,:), Forecast, ForecastCov)


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
base_periods = 9;
T = 24; 

% Historical data, fill NaN data with last value
Xorig = table2array(dataUS(:,2:end));
Xorig(end,1) = Xorig(end-1,1);          % assume zero growth IP in 3/2016

% Baseline 2016
Xbase = zeros(base_periods,size(Xorig,2));

% X(:,1) Industrial production 
% growth by 2 percent per annum during 2016
IndProdGrowthAnnum = 0.02;
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
shock = 1;
Xbase(end,3) = shock + Xorig(end,3);
step = (Xbase(end,3) - Xorig(end,3))/base_periods;
Xbase(:,3) = cumsum(repmat(step,base_periods,1))+Xorig(end,3);

% X(:,4) US 10-year yield
% 2016 200 bps hike
shock = 2;
step = shock/base_periods;
Xbase(:,4) = cumsum(repmat(step,base_periods,1))+Xorig(end,4);

% X(:,5) VIX
% 2016 Up by 1 standard deviation
stdVIX = std(dataUS.VIX);
shock = stdVIX;
step = shock/base_periods;
Xbase(:,5) = cumsum(repmat(step,base_periods,1))+Xorig(end,5);

% X(:,6) Dow Jones
% 2016 Down by 1 standard deviation, growth rate
stdDow = -std(dY_us12(:,6));
step = (1+stdDow)^(1/base_periods);
Xbase(:,6) = cumprod(repmat(step,base_periods,1))*Xorig(end,6);

% Extend historical series by 2016 baseline
% calculate 1-month growth rates
Xextend = [Xorig; Xbase];
dXextend = diff(log(Xextend));

% 24-month forecast for the remaining values
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
stdshock = 1;  % standard deviation shock, all except IGAE
stdIGAE  = 2;  % standard deviation shock, IGAE

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
idx = find(sum(idx,1)>=20);


%% Loose end - path of US variables

% 12-month path of US variables

dY94 = dY_us01(33:44,:);    % select US initial baseline 10/94 - 10/95
base_len = size(dY94,1);    % length, initial baseline
dY94 = [dY_us01; dY94];     % add baseline scenario to previous
lagvector = [0 1 2 3 4];
X = lagmatrix(dY94, lagvector);

idxHistPlusBase = [idxMXnotNaN;(1:1:base_len)' + idxMXnotNaN(end)];
X=X(idxHistPlusBase,:);

% Cell structure, design matrix
Xexpand = kron(X,eye(nMX));
Xcell = mat2cell(Xexpand, nMX*ones(size(X,1),1));
nX = size(Xexpand,2);




%% Section IV - create shocks

% Create scenario tables


months = 1:1:12;
years = 2016:1:2020;
nmonths = 12;
nyears = 5;

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
scenario_array = zeros(nmonths*nyears,size(dataUS,2)-1);

% create data table

dataScenario = table(scenario_dates,scenario_array(:,1), ...
    scenario_array(:,2),scenario_array(:,3),scenario_array(:,4), ...
    scenario_array(:,5),scenario_array(:,6), ...
    'VariableNames',{'Date','IndProd','WTI','yld3mo','yld10yr','VIX','DJ'});

dataScenario = dataScenario(3:end,:);  % eliminate first three months of 2016


%% Loose ends
auxspec = vgxvarx(mdlUS, Y(5:end,:),[],Y(1:4,:));
T=36;
[Forecast, ForecastCov] = vgxpred(VARmdlUS,T,[], Y(end-3:end,:));
vgxplot(VARmdlUS, Y(end-12:end,:), Forecast, ForecastCov)














