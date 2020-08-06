close all;
clear all;
%%---VARIABLES---
interval = linspace(180,1000,1000);
sigma=[1.8000e+2; 5.2000e+2; 7.2000e+2; 7.5000e+2; 8.0000e+2; 1.0000e+3];

e1=[5.0000e-04;1.3000e-03;2.0000e-03;4.5000e-03;6.0000e-03;8.5000e-03];
e2=[1.0000e-03;1.8000e-03;2.5000e-03;2.6197e-03;2.8276e-03;3.7655e-03];
e3=[1.2000e-03;2.0000e-03;3.7000e-03;3.8197e-03;4.0276e-03;5.9655e-03];

sigmap=735;
lengthm=10;
n=length(sigma);
xinterval=linspace(7.2000e+2,7.5000e+2,1000);
%---------------------------------------------------

%---STRESS EXPERIMENT FOR EPSILON 1
[p1,s1,mu1]=polyfit(sigma,e1,n-1);
p1y=polyval(p1,interval,s1,mu1);
s1y=spline(sigma,e1,interval); 
[b01,b11]=linearapprox(sigma,e1,n);
y1=lineareval(b01,b11,interval);

%EVALUATE SIGMAP
e1p=polyval(p1,sigmap,s1,mu1);
eps1=spline(sigma,e1,sigmap);
yep1=lineareval(b01,b11,sigmap);

%MIN AND MAX VALUE IN THE INTERVAL xi<sigmap<xi+1
pval1=polyval(p1,xinterval,s1,mu1);
sval1=spline(sigma,e1,xinterval);
lval1=lineareval(b01,b11,xinterval);
minP1 = min(pval1);
maxP1 = max(pval1);
minS1 = min(sval1);
maxS1 = max(sval1);
minL1 = min(lval1);
maxL1 = max(lval1);

%PLOT THE INTERPOLATION GRAPHS
figure
subplot(3,1,1);
plot(sigma,e1,'*r')
hold on
plot(interval,p1y,'-g');
plot(interval,s1y,'-b');
plot(interval, y1, '-k');
hold off
legend('\epsilon1','Pol','Spline','LinearApprox','Location','NorthWest');
title('Subplot 1: STRESS EXPERIMENT FOR \epsilon1');
xlabel('\sigma');
ylabel('f(\sigma)');
ax = gca;
ax.FontSize = 11;
%-------------------------------------------------

%---STRESS EXPERIMENT FOR EPSILON 2
[p2,s2,mu2]=polyfit(sigma,e2,2);
p2y=polyval(p2,interval,s2,mu2);
s2y=spline(sigma,e2,interval);
[b02,b12]=linearapprox(sigma,e2,n);
y2=lineareval(b02,b12,interval);

%EVALUATE SIGMA P
e2p=polyval(p2,sigmap,s2,mu2);
eps2=spline(sigma,e2,sigmap);
yep2=lineareval(b02,b12,sigmap);

%MIN AND MAX VALUE IN THE INTERVAL xi<sigmap<xi+1
pval2=polyval(p2,xinterval,s2,mu2);
sval2=spline(sigma,e2,xinterval);
lval2=lineareval(b02,b12,xinterval);
minP2 = min(pval2);
maxP2 = max(pval2);
minS2 = min(sval2);
maxS2 = max(sval2);
minL2 = min(lval2);
maxL2 = max(lval2);

%PLOT THE INTERPOLATION GRAPHS
subplot(3,1,2);
plot(sigma,e2,'*r')
hold on
plot(interval,p2y,'-g');
plot(interval,s2y,'-b');
plot(interval,y2,'-k');
hold off
legend('\epsilon2','Pol','Spline','LinearApprox','Location','NorthWest');
title('Subplot 2: STRESS EXPERIMENT FOR \epsilon2')
xlabel('\sigma');
ylabel('f(\sigma)');
ax = gca;
ax.FontSize = 11;
%-------------------------------------------

%---STRESS EXPERIMENT FOR EPSILON 3
[p3,s3,mu3]=polyfit(sigma,e3,n-1);
p3y=polyval(p3,interval,s3,mu3);
s3y=spline(sigma,e3,interval);
[b03,b13]=linearapprox(sigma,e3,n);
y3=lineareval(b03,b13,interval);
%EVALUATE SIGMAP
e3p=polyval(p3,sigmap,s3,mu3);
eps3=spline(sigma,e3,sigmap);
yep3=lineareval(b03,b13,sigmap);

%MIN AND MAX VALUE IN THE INTERVAL xi<sigmap<xi+1
pval3=polyval(p3,xinterval,s3,mu3);
sval3=spline(sigma,e3,xinterval);
lval3=lineareval(b03,b13,xinterval);
minP3 = min(pval3);
maxP3 = max(pval3);
minS3 = min(sval3);
maxS3 = max(sval3);
minL3 = min(lval3);
maxL3 = max(lval3);

%PLOT THE INTERPOLATION GRAPHS
subplot(3,1,3);
plot(sigma,e3,'*r')
hold on
plot(interval,p3y,'-g');
plot(interval,s3y,'-b');
plot(interval,y3,'-k');
hold off
legend('\epsilon3','Pol','Spline','LinearApprox','Location','NorthWest');
title('Subplot 3: STRESS EXPERIMENT FOR \epsilon3')
xlabel('\sigma');
ylabel('f(\sigma)');
ax = gca;
ax.FontSize = 11;

%linear approximation
function [b0,b1]=linearapprox(xi,yi,n)
    s1=sum(xi);
    xsq=xi.^2;
    s2=sum(xsq);
    v1=sum(yi);
    xy=xi.*yi;
    v2=sum(xy);
    D=((n+1)*s2)-(s1^2);
    b0=(1/D)*((v1*s2)-(v2*s1));
    b1=(1/D)*(((n+1)*v2)-(v1*s1));
end
function px=lineareval(b0,b1,x)
    px=b0+(b1*x);
end

