function [X,Y]=find_XY(y,p)
% Change VAR(p) model to VAR(1)
% Inputs
% y : T * K data
% p : scalar, number of lags
% output
% X :  T-p * kp matrix , stack of lags
% X :  T-p * k  matrix , dependent variable

[T,k]=size(y);
X=zeros(T-p,k*p);                   %initialize variable
for i=1:p
    X(:,i*k-k+1:i*k)=y(p-i+1:T-i,:);
end
Y=y(p+1:end,:);