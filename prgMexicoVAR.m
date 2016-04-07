%% prgMexicoVAR
%
%  Author: Jorge A. Chan-Lau
%  Current version: 04/06/2016
%  Previous version: 
%  
%  Requires Matlab Econometrics Toolbox (different than Le Sage toolbox)
%  Tested in Matlab 2015a, may run in Matlab 2014b
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

%% Section II: VAR specification 

% US model
% --------
nUS = size(dY_us01,2);
nAR = 4;
cons = true;
Y = dY_us01;
idxUSnotNaN = find(~isnan(sum(Y,2)));   % used only valid observations
Y = Y(idxUSnotNaN,:);


% set the model
mdlUS = vgxset('n',nUS,'nAR',nAR,'Constant',cons);

[VARmdlUS, VARmdlUSStdErrors, logLmdlUS, W] = ...
    vgxvarx(mdlUS,Y(5:end,:),[],Y(1:4,:));

auxspec = vgxvarx(mdlUS, Y(5:end,:),[],Y(1:4,:));
T=36;
[Forecast, ForecastCov] = vgxpred(VARmdlUS,T,[], Y(end-3:end,:));
vgxplot(VARmdlUS, Y(end-12:end,:), Forecast, ForecastCov)

% Mexican model
% -------------

nMX = size(dY_mx01,2);
nAR = 4;
cons = true;
Y = dY_mx01;
idxMXnotNaN = find(~isnan(sum(Y,2)));   % used only valid observations
Y = Y(idxMXnotNaN,:);
X = dY_us01(idxMXnotNaN,:);

Z = [X Y];  % expanded matrix

% Model setup
% Assumes that US variables do not depend on MX variables

% 
% A generic AR(p) matrix is:
%
%    AR(US-US)(6x6) |       0 (6x7)
%   --------------------------------
%    AR(US-MX)(7x6) | AR(MX-MX)(7x7)
%

% Create logical AR(p) matrices

logical_AR = true(nUS+nMX);
logical_AR(1:nUS,(nUS+1):(nUS+nMX)) = false;
num_AR = double(logical_AR);
ARsolve = repmat({logical_AR},nAR,1);
AR = repmat({num_AR},nAR,1);

% set model with vgxset

mdlMX = vgxset('n',nUS+nMX,'nAR',nAR,'Constant',cons, ...
    'AR',repmat({num_AR},nAR,1), ...
    'ARsolve',repmat({logical_AR},nAR,1));

[VARMX, VARMXStdErrors, logLMX, WMX] = ...
    vgxvarx(mdlMX,Z(5:end,:),[],Z(1:4,:));

auxspec = vgxvarx(mdlUS, Y(5:end,:),[],Y(1:4,:));
T=36;
[Forecast, ForecastCov] = vgxpred(VARmdlUS,T,[], Y(end-3:end,:));
vgxplot(VARmdlUS, Y(end-12:end,:), Forecast, ForecastCov)


%% Section III shock specification




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
















