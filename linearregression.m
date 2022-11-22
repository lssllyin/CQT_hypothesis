function [a1,b1,r2,p]=linearregression(XX,YY,C,lspec)
%%% XX is the independent variables, YY is the dependent variables,C is the
%%% color of regression line
%%% example: [a,b,r,p]=linearregression(1./CO2_20(r5),iso_20(r5),'--r');
x=XX;
y=YY;
color=C;
n = length(x);
x1=[ones(n,1) x];
[b,bint,r,rint,stats]=regress(y,x1);
% RMSE=sqrt(nanmean(r.^2));
r2=stats(1);
p=stats(3);
a1=b(1);%intercpt
b1=b(2);% slope
% v=axis;
X=min(x):(max(x)-min(x))/n:max(x);
Y=b(1)+b(2).*X;
plot(X,Y,'color',color,'LineWidth',2,'Linestyle',lspec);
end