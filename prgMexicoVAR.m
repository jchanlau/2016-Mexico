%% prgMexicoVAR
%
%  Author: Jorge A. Chan-Lau
%  Current version: 04/06/2016
%  Previous version: 
%  
%  Requires Matlab Econometrics Toolbox (different than Le Sage toolbox)
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
%
%  US variables
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
%
%  Mexico variables







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

[VARmdlUS, VARmdlUSStdErrors, logLmdlUS, W] = vgxvarx(mdlUS,dY_us01);






%% Section III shock specification









