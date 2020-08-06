close all
clear all

%Script to read the COVID-19 Table for Marche region
TT=readtable('covid19_hospitalized_italy.csv'); % TT is a table structure

y=TT.Marche;     % y is a vector
t=datenum(TT.date); % transform dates in a numerical format
%---END SCRIPT---%

%---VARIABLES---
yN = length(y);
tN = length(t);
Ne = 6;
Np = 20;
tint=t(1:(tN-Ne));
yint=y(1:(yN-Ne));
text=t((tN-Ne)+1:tN);
yext=y((tN-Ne)+1:yN);

%%
%---Interpolation using POLYNOMIAL function---%
enint=zeros(Np,1);
enext=zeros(Np,1);
for n=1:Np
    x = linspace(t(1),t(tN-Ne),n+1);
    x=round(x,0);
    r = linspace(1,t(tN-Ne)-t(1),n+1);
    r=round(r,0);
    for i=1:n+1
        h(i) = y(r(i));
    end
    [Px,s,mu] = polyfit(x,h,n);
    Pn = polyval(Px,t,s,mu);
    for i=3:(tN-Ne)
        Perr(i) = polyval(Px,tint(i),s,mu);
        enint(n)=enint(n)+norm((Perr(i)-yint(i)))/norm(yint(i));
    end
    for i=1:Ne
        Perr2(i) = polyval(Px,text(i),s,mu);
        enext(n)=enext(n)+norm((Perr2(i)-yext(i)))/norm(yext(i));
    end
    minPn(n)=min(Pn);
    maxPn(n)=max(Pn);
end

%FIND THE DEGREE OF THE POLYNOMIAL WITH MINIMIZE THE ERROR
[minInt,dint]=min(enint);
[minExt,dext]=min(enext);

%COMPUTE THE POLYNOMIAL
x1 = linspace(t(1),t(tN-Ne),dint+1);
x2 = linspace(t(1),t(tN-Ne),dext+1);
x1=round(x1,0);
x2=round(x2,0);
r1 = linspace(1,t(tN-Ne)-t(1),dint+1);
r1=round(r1,0);
r2 = linspace(1,t(tN-Ne)-t(1),dext+1);
r2=round(r2,0);
for i=1:dint+1
    h1(i) = y(r1(i));
end
for i=1:dext+1
    h2(i) = y(r2(i));
end
[Px1,s1,mu1] = polyfit(x1,h1,dint);
[Px2,s2,mu2] = polyfit(x2,h2,dext);
P1 = polyval(Px1,tint,s1,mu1);
P2 = polyval(Px2,t,s2,mu2);


%Plot the data
figure
plot(t,y,'*r')
hold on
plot(tint,P1,'-g')
plot(t,P2,'-b')
datetick('x','mmm-dd')
legend('Hosp','P_{int}','P_{ext}');
title('Marche (Best Polynomial)')
ylabel('Hospitalized individuals')
xlabel('Date');

ylim([-100 1300]);
%--------------------------------------------------------------
%---SPLINE---%
clear h;
eensn=zeros(Np,1);
evnsn=zeros(Np,1);
for n=1:Np
    x = linspace(t(1),t(tN-Ne),n+1);
    x=round(x,0);
    r = linspace(1,t(tN-Ne)-t(1),n+1);
    r=round(r,0);
    for i=1:n+1
        h(i) = y(r(i));
    end
    pp1 = spline(x,h);
    Sn = ppval(pp1,t);
    for i=3:(tN-Ne)
        Serr(i) = ppval(pp1,tint(i));
        eensn(n)=eensn(n)+norm((Serr(i)-yint(i)))/norm(yint(i));
    end
    for i=1:Ne
        Serr2(i) = ppval(pp1,text(i));
        evnsn(n)=evnsn(n)+norm((Serr2(i)-yext(i)))/norm(yext(i));
    end
    minSn(n)=min(Sn);
    maxSn(n)=max(Sn);
end

[minInt,dint]=min(enint);
[minExt,dext]=min(enext);

%COMPUTE THE POLYNOMIAL
x1 = linspace(t(1),t(tN-Ne),dint+1);
x2 = linspace(t(1),t(tN-Ne),dext+1);
x1=round(x1,0);
x2=round(x2,0);
r1 = linspace(1,t(tN-Ne)-t(1),dint+1);
r1=round(r1,0);
r2 = linspace(1,t(tN-Ne)-t(1),dext+1);
r2=round(r2,0);
for i=1:dint+1
    h1(i) = y(r1(i));
end
for i=1:dext+1
    h2(i) = y(r2(i));
end
ppx1 = spline(x1,h1);
ppx2= spline(x2,h2);
S1 = ppval(ppx1,tint);
S2 = ppval(ppx2,t);

%---PLOT THE RESULTS---
figure
plot(t,y,'or')
hold on
plot(tint,S1,'-g');
plot(t,S2,'-k');
datetick('x','mmm-dd')

legend('Hosp','S_{int}','S_{ext}');
title('Marche (Best Spline)')
ylabel('Hospitalized individuals')
xlabel('t');
ylim([-100 1300]);
ax = gca;
ax.FontSize = 11;
%-----------------------------------------------------------
%%
%%---PART2---
%%--ASSIGNMENT 1
[pl,s,mu] = polyfit(t(3:tN-Ne),log(y(3:yN-Ne)),1);
alfa=pl(1);
Plin=polyval(pl,t,s,mu);

%---COMPUTE THE RELATIVE ERRORS---
for i=3:(tN-Ne)
    Perr(i)=polyval(pl,t(i),s,mu);
    reint(i)=norm(Perr(i)-log(y(i)))/norm(log(y(i)));
end
for i=1:Ne
    Perr2(i)=polyval(pl,text(i),s,mu);
    reext(i)=norm(Perr2(i)-log(yext(i)))/norm(log(yext(i)));
end

%--PLOT THE RESUTS----

figure
plot(tint(3:(tN-Ne)),log((yint(3:(tN-Ne)))),'or')
hold on
plot(t,Plin,'-g')
datetick('x','mmm-dd')

legend('Hosp','Linear Regression');
title('Marche (Linear Regression)')
ylabel('log(Hospitalized individuals)')
xlabel('t');

ax = gca;
ax.FontSize = 11;
%%
%%---ASSIGNMENT 2---
i=3;
reint2 = zeros(8,1);
reext2 = zeros(Ne,1);
k=1;
m=1;
NW= fix((tN-Ne)/6)*6;
while i<=NW
    d6=t(i:i+5); %6 days window
    h6=y(i:i+5); %hospitalized in 6 days
    [px6,s6,mu6]=polyfit(d6,log(h6),1);
    alpha6(k)=px6(1);
    if i==3
        pxw1=px6;
        sw1=s6;
        muw1= mu6;
    end
    for j=1:6
        p6(j)=polyval(px6,d6(j),s6,mu6);
        reint2(k)=reint2(j)+norm(p6(j)-h6(j))/norm(h6(j)); %relative errors
        v(m)=p6(j); 
        m=m+1;
    end
    i=i+6;
    k=k+1;
end

%LAST WINDOW OF 2 DAYS
d6=t((NW+3):tN-Ne); %last window
h6=y((NW+3):tN-Ne); %last hospitalized
l=length(p6);
[px6,s6,mu6]=polyfit(d6,log(h6),1);
alpha6(k)=px6(1);
for j=1:(tN-Ne-(NW+2))
    p6(j)=polyval(px6,d6(j),s6,mu6);
    reint2(k)=reint2(j)+norm(p6(j)-h6(j))/norm(h6(j)); %relative errors
    v(m)=p6(j); 
    m=m+1;
end
k=k+1;
%------------------------------------------------------------------
%---RELATIVE ERROR USING THE LAS WINDOW POLYNOMIAL FOR THE EXTERNAL 6 POINTS
d6=t(tN-Ne+1:tN); %external 6 points
h6=y(tN-Ne+1:yN); %last hospitalized
for j=1:Ne
    p6(j)=polyval(px6,d6(j),s6,mu6);
    Perr2(j)=polyval(px6,text(j));
    reext2(k)=reext2(j)+norm(Perr2(j)-yext(j))/norm(yext(j));
    v(m)=p6(j); 
    m=m+1;
end

figure
plot(t,log(y),'or')
hold on
plot(t(3:8),v(1:6),'-k')
plot(t(9:14),v(7:12),'-k')
plot(t(15:20),v(13:18),'-k')
plot(t(21:26),v(19:24),'-k')
plot(t(27:32),v(25:30),'-k')
plot(t(33:38),v(31:36),'-k')
plot(t(39:44),v(37:42),'-k')
plot(t(45:46),v(43:44),'-k')
plot(t(47:52),v(45:50),'-k')
datetick('x','mmm-dd')

legend('Hosp','Linear Regression');
title('Marche (6 Days Window)')
ylabel('log(Hospitalized individuals)')
xlabel('t');
ax = gca;
ax.FontSize = 11;
%%
%%---ASSIGNMENT 3
Pw1 = polyval(pxw1,t,sw1,muw1);
figure
plot(t,log(y),'*r')
hold on
plot(t,Pw1,'-k')
yline(log(1525271),'-b')
datetick('x','mmm-dd')
legend('Hosp','Linear Regression','Marche population');
title('Marche (Linear Regression with \alpha_1)')
ylabel('log(Hospitalized individuals)')
xlabel('t');
ax = gca;
ax.FontSize = 11;