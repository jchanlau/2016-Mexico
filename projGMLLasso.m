function [ mdl ] = projGMLLasso(Y,X,nlag)
% Lasso estimation 
%   Y       NPL, to be transformed into logits
%   X       matrix of covariates
%   nlag    desired number of lags

% create lag matrix

[B,FitInfo]=lassoglm(Zlag,Y,'binomial','NumLambda',25,'CV',10);

%[B,FitInfo]=lassoglm(Zlag,Y,'normal','NumLambda',25,'CV',10);

lassoPlot(B,FitInfo,'PlotType','CV');
lassoPlot(B,FitInfo,'PlotType','Lambda','Xscale','log');

idx = FitInfo.Index1SE;
B0=B(:,idx);
cnst = FitInfo.Intercept(idx);
B1 = [cnst;B0];

prds = glmval(B1,Zlag,'logit');
plot([prds Y]);
prds1 = glmval(B1,
end

