% PDprojection.m
% Author: Jorge A. Chan-Lau
%
%
% Machine learning selection, predictive models
% Data input: 
%      CNVBData.xlsx                bank data
%      MacroPredictors.xlsx         average macro paths
%           Sheets: Adverse, Severe
% -------------------------------------------------------------------------            

clear all; clc;
bank_file = 'CNVBData.xlsx';
bank_name = {
    'BBVA Bancomer';
    'Banamex';
    'Santander';
    'Banorte';
    'HSBC';
    'Inbursa';
    'ScotiaBank';
    };

bank_number = size(bank_name,1);

%% READ COVARIATES

macrofile = 'MacroPredictors.xlsx';
macrosheet = 'Adverse';
[Z,B] = xlsread(macrofile, macrosheet);

% find indices initial date and final date
init_date = '1/31/2001 ';
end_date = '12/31/2015';
init_dateidx = strmatch(init_date,B(2:end,:));
end_dateidx = strmatch(end_date,B(2:end,:));

% create indices for training and testing datasets
train_idx = init_dateidx:end_dateidx';
test_idx =(end_dateidx+1):size(Z,1);

%Covariates
dZ = [repmat(nan,1,size(Z,2));diff(log(Z))];
dZtrain = dZ(train_idx,:);
dZtest = dZ(test_idx,:);
dZtot = dZ([train_idx test_idx],:);
% Z-scores of Z variables, in levels
idx_notNaN = max(find(isnan(sum(Z,2))))+1;
Znorm = zscore(Z(idx_notNaN:end,:));
Znorm = [repmat(NaN,idx_notNaN-1,size(Znorm,2)); Znorm];



%% ESTIMATE BANK MODELS



for i = 1:bank_number,
    disp(bank_name(i));
    
    % Read Bank data    
    
    [bank_data, B] = xlsread(bank_file,char(bank_name(i)));
    
    %----------------------------------------------------------------------
    % Bank total assets, in mn local currency
    bank_assets = bank_data(1,:);
    
    %----------------------------------------------------------------------
    % Loans ( in mn local currency)
    
    loan_total = bank_data(11,:);               % model estimation
    % Corporate and gov/loans
    loan_corp_total  = bank_data(12,:);         % model estimation
    loan_corp_nonfin = bank_data(13,:);         % model estimation
    loan_corp_fin = bank_data(20,:);            % model estimation
    loan_corp_gov = bank_data(24,:);            % model estimation
    loan_corp_gov_fed = bank_data(25,:);
    loan_corp_gov_mun = bank_data(26,:);
    loan_corp_gov_noguarantee = bank_data(27,:);
    % Households
    loan_hh_total = bank_data(29,:);            % model estimation
    loan_hh_creditcard = bank_data(30,:);       
    loah_hh_nonrecurring = bank_data(31,:);
    loan_hh_durable = bank_data(34,:);
    % Mortgages
    loan_mtg_total = bank_data(39,:);           % model estimation
    loan_mtg_resid = bank_data(40,:);           % model estimation
    loan_mtg_lowincome = bank_data(41,:);       % model estimation
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Investments (in mn local currency)
    
    assets_total = bank_data(7,:);
    assets_trade = bank_data(8,:);
    assets_sale  = bank_data(9,:);
    assets_book = assets_total - assets_trade - assets_sale;
    
    %----------------------------------------------------------------------
    % Non-performing loans (in percent)
    
    npl_total = bank_data(195,:);               % model estimation
    % Corporate and gov/loans
    npl_corp_total = bank_data(196,:);          % model estimation
    npl_corp_nonfin = bank_data(197,:);         % model estimation
    npl_corp_fin = bank_data(198,:);            % model estimation
    npl_corp_gov = bank_data(199,:);            % model estimation
    % Households    
    npl_hh_total = bank_data(202,:);            % model estimation
    npl_hh_creditcard = bank_data(203,:);       
    npl_hh_nonrecurring = bank_data(204,:);
    npl_hh_durable = bank_data(207,:);
    % Mortgages
    npl_mtg_total = bank_data(212,:);           % model estimation
    npl_mtg_resid = bank_data(213,:);           % model estimation
    npl_mtg_lowincome = bank_data(214,:);       % model estimation
    % ---------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Income components (in mn local currency)
       
    net_interest_income = bank_data(144,:);
    h = reshape(net_interest_income(2:end),12,15);
    h(2:end,:) = h(2:end,:) - h(1:end-1,:);
    h=reshape(h,1,180);
    net_interest_income(2:end) = h;
    
    net_income_beftax = bank_data(185,:);
    h = reshape(net_income_beftax(2:end),12,15);
    h(2:end,:) = h(2:end,:) - h(1:end-1,:);
    h=reshape(h,1,180);
    net_income_beftax(2:end) = h;    
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Group only relevant variables
    % Eliminate first observation,December 2000
    %
    
    Loans = [ ...
        loan_total;             % 1: total
        loan_corp_total;        % 2: commercial, total
        loan_corp_nonfin;       % 3: commercial, nonfinancials
        loan_corp_fin;          % 4: commercial, financials
        loan_corp_gov;          % 5: commercial, government
        loan_hh_total;          % 6: households, total
        loan_mtg_total;         % 7: mortgages, total
        loan_mtg_resid;         % 8: mortgages, high income/middle class
        loan_mtg_lowincome;     % 9: mortgages, low income
        ];
    Loans2assets = bsxfun(@ldivide,bank_assets, Loans);
    Loans = Loans(:,2:end);
    Loans2assets = Loans2assets(:,2:end);
    
    Income = [ ...
        net_interest_income;
        net_income_beftax; 
        ];
    Income2assets = bsxfun(@ldivide, bank_assets, Income);
    Income = Income(:,2:end);
    Income2assets = Income2assets(:,2:end);
        
    NPL = [ ...
        npl_total; 
        npl_corp_total;
        npl_corp_nonfin;
        npl_corp_fin;
        npl_corp_gov;
        npl_hh_total;
        npl_mtg_total;
        npl_mtg_resid;
        npl_mtg_lowincome;
        ];
    dNPL = diff(log(NPL/100)');    
    NPL=NPL(:,2:end)'/100;
    
    %----------------------------------------------------------------------    
    % REGULARIZED GLM ESTIMATIONS
    
     
    % range for estimation    
    sample_all = 1:1:size(NPL,1);
    sample_third= round(size(NPL,1)/3):1:size(NPL,1);
    sample_half = round(size(NPL,1)/2):1:size(NPL,1);
    my_sample = sample_all;
  
    % NPL predictions
    n_loans=size(NPL,2);
    loan_idx = [1 2 3 6 7];
    predicted_NPL = zeros(size(Znorm,1),length(loan_idx))+ NaN;
    
    glmcell = cell(length(loan_idx),1);
    glmfit  = cell(length(loan_idx),1);

    for j = 1:length(loan_idx),
        disp(j);
        var_loan = loan_idx(j);        
        nlag = 4;
        X = Znorm(:,7:end);
        X = lagmatrix(X,0:nlag);
        cnst = ones(size(X,1),1);
        X = [cnst,X];        
        rng('default');
        Y = NPL(:,var_loan);
        XX = X(train_idx,:);
        new_idx = max(find(isnan(Y)));
        if ~isempty(new_idx),
            Y = Y(new_idx+1:end,:);
            XX = XX(new_idx+1:end,:);            
        end,
        if length(my_sample)<size(XX,1) 
            Y = Y(my_sample,:);
            XX = XX(my_sample,:);           
        end, 
        [B,FitInfo]=lassoglm(XX,Y,'binomial','NumLambda',25,'CV',10);   
%        lassoPlot(B,FitInfo,'PlotType','CV');
%        lassoPlot(B,FitInfo,'PlotType','Lambda','Xscale','log');
        idx = FitInfo.Index1SE;
        B0=B(:,idx);
        intcpt = FitInfo.Intercept(idx);
        B1 = [intcpt;B0];
        prdhist = glmval(B1,XX,'logit');
        %plot([prdhist Y]);
        Xsim = X(test_idx,:);
        prdsim = glmval(B1,Xsim,'logit');
        prds = [prdhist; prdsim];
        predicted_NPL(end-(size(prds,1))+1:end,j) = prds;        
        glmcell{j}=B1;
        glmfit{j} = FitInfo;
    end,
    
    bank(i).name = bank_name(1);
    bank(i).predNPL = predicted_NPL;
    bank(i).NPL = NPL(:,loan_idx);
    Loans = Loans';
    bank(i).Loans = Loans(:,loan_idx)';
    bank(i).Fit = glmfit;
    bank(i).glmcoef = glmcell;
    
end,
Loans = Loans';
Loans = Loans(:,loan_idx);

save('PDprojection.mat',...
    'bank',...
    'test_idx','train_idx',...
    'Loans'...
    );


    


%% Loose ends

% Net income prediction
    Y=Income2assets(2,:)';
    X = dZ(:,7:end);
    X = lagmatrix(X,0:nlag);
    cnst = ones(size(X,1),1);
    X = [cnst,X];        
    rng('default');
    XX = X(train_idx,:);
    Y = Y(my_sample,:);
    XX = XX(my_sample,:);
    [B,FitInfo]=lasso(XX,Y,'NumLambda',25,'CV',10);   
    lassoPlot(B,FitInfo,'PlotType','CV');    
    lassoPlot(B,FitInfo,'PlotType','Lambda','Xscale','log');
    idx = FitInfo.Index1SE;
    B0=B(:,idx);
    intcpt = FitInfo.Intercept(idx);
    B1 = [intcpt;B0];
    prdhist = glmval(B1,XX,'normal');
        %plot([prdhist Y]);
        Xsim = X(test_idx,:);
        prdsim = glmval(B1,Xsim,'logit');
        prds = [prdhist; prdsim];

