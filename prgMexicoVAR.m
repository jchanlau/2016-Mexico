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
%  17.   IGAE
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

%% Section II: VAR specification 

% US model
nUS = size(dY_us01,2);
nAR = 4;
cons = true;

% set the model
mdlUS = vgxset('n',nUS,'nAR',nAR,'Constant',cons);

[VARmdlUS, VARmdlUSStdErrors, logLmdlUS, W] = ...
    vgxvarx(mdlUS,dY_us01(5:end,:),[],dY_us01(1:4,:));

T=10;
[Forecast, ForecastCov] = vgxpred(VARmdlUS,T,[], dY_us01(end-9:end,:));
vgxplot(VARmdlUS, dY_us01(end-9:end,:), Forecast, ForecastCov)


%% Section III shock specification









