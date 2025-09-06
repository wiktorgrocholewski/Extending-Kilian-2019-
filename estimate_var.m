function [A_hat,Sigma_hat_u,u_hat] = estimate_var(y,p)
%Estimates Var(p) model
% Inputs
% y : T * K data
% p : scalar, number of lags

%output 
%A_hat : k*kp estimate of Var coef


[X,Y]= find_XY(y,p);             % change the model to VAR(1)

[T,k]=size(Y);
%A_hat1= kron(eye(k),X)\Y(:);     %Apply OLS
%A_hat=reshape(A_hat1,[],k)';
A_hat = (X' * X) \ (X' * Y);
A_hat = A_hat';


u_hat=Y-X*A_hat';                   %Find resduals
% Find the estimated error covariance matrix 
Sigma_hat_u=u_hat'*u_hat/T;           


