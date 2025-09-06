clear all
close all
clc


% This file is meant for SVAR identification using heteroskedasticity
% We used it as the main file for our assignment


% First I load the data and select only the variables of interest

data = readtable('ourdata.csv');
data = data(:,1:6);
dataT = table2array(data);
[len, numvars] = size(dataT);

addpath '/Users/wiktor/School/Tilburg/Adv Macrometrics/Assignment Hetero Code/Mbreak_matlab'


% Estimation sample
smplStart = '1974M01';
smplEnd   = '2017M12'; 


varNames_paper = {'Real oil price','World oil production','World oil inventories','Global real activity','U.S. industrial production','U.S. CPI'};
varNames_paperVD = {'Real oil price','Oil production','Oil inventories','World IP','U.S. IP','U.S. CPI'};

% I used 24 lags, to get comparable output with Kilian (2023)

lags = 24;

dep_vars = dataT(lags+1:end,:);
regressors = zeros(len-lags, numvars*lags);
for c = 1:lags
    lagged = dataT(lags+1-c:end-c,:);
    regressors(:,6*c-5:6*c) = lagged;
end
regressors = [ones(len-lags,1)';regressors']';



% estimating reduced form VAR

[A_hat,Sigma_hat_u,u_hat] = estimate_var(dataT,lags);  %Estimating the Model

% Dividing before and after 208

u_hat1 = u_hat(1:206,:);
T1 = rows(u_hat1);
u_hat2 = u_hat(207:end,:);
T2 = rows(u_hat2);

Sigma_hat_u1=u_hat1'*u_hat1/T1;  
Sigma_hat_u2=u_hat2'*u_hat2/T2;  

[D0_hat, lambda] = decompose_hetero(Sigma_hat_u1, Sigma_hat_u2);

trans = ones(6,1)*[1 1 1 1 1 -1];

%We do a sign flip for 6th shock

D0_hat = D0_hat.*trans;

H=50;
Theta_hat  = irf(A_hat,D0_hat,H);  



% compute I.I.D. bootstrap CI

B = 399;
[X,Y]=find_XY(dataT,lags);
[A_hat,sigma_hat_u,u_hat] = estimate_var(dataT,lags);




% Important note: In the end for our paper we used a t-percentile
% bootstrap, however with 24 lags it takes almost 3 hours to run on my
% computer, so in this verion of the code I use naive bootstrap. To change
% it at the end of the inputs for the fixedreg_boot one needs to change
% 'naive' for 't-percentile'



% IRF fixed regressor bootstrap with naive or t-percentile CI

[CI05,CI16,CI84,CI95] = fixedreg_boot(X,A_hat,u_hat1,u_hat2,Sigma_hat_u1,Sigma_hat_u2,lags,Theta_hat,D0_hat,B,H,'naive');

% Scaling so that it matches Kilian (2023)
scale_param = 10 / Theta_hat(1,6,1);

Theta_hat = Theta_hat * scale_param;
CI05 = CI05 * scale_param;
CI95 = CI95 * scale_param;
CI16 = CI16 * scale_param;
CI84 = CI84 * scale_param;


% Last argument is about the structural shock that we want to plot, so in our
% case after analysis we decided that 6th shock is most likely the Oil
% Price shock that we were trying to estimate.
plot_irf(Theta_hat,CI05,CI16,CI84,CI95,6)
% sgtitle('IRF to Oil Price Shock with naive or t-percentile CI', 'FontSize', 30, 'FontWeight', 'bold');









