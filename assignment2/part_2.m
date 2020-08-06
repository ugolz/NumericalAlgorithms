close all
clear all


TT=readtable('covid19_regional_data.xlsx','Sheet','Dead'); % TT is a table structure
data_dead=TT.Marche;     % y is a vector
t=datenum(TT.date); % transform dates in a numerical format
clear TT

TT=readtable('covid19_regional_data.xlsx','Sheet','NewPositive'); % TT is a table structure
data_positive=TT.Marche;     % y is a vector
clear TT

TT=readtable('covid19_regional_data.xlsx','Sheet','Recovered'); % TT is a table structure
data_recovered=TT.Marche;     % y is a vector
clear TT
% Parameters for the Cauchy problem
t0= datenum('2/1/2020');  % initial time
T=  t(length(t));%t0+100;               % final time

%% parametri SIR

par(1)=1525271; %Np %ESTIMATION OF MARCHE POPULATION
par(2)=2;  %beta0
par(3)=1/5;  %delta;
par(4)=0.4; %sigma;
par(5)=1/7; %gammaA;
par(6)=1/14; %gammaR;
par(7)=1/21; %alpha;
par(8)= 0.2;  %betaA;

%% initial condition
y0=zeros(6,1);

y0(2)=1;   %E   
y0(1)=par(1)-y0(2); %S       
y0(3)=0;   %A
y0(4)=0;   %I
y0(5)=0;   %R
y0(6)=0;   %D

varname={'S','E','A','I','R','D'}; % names for figures

% solution using ode45
h=1;
time_val= t0:h:T;

% MATLAB SOLVER
options=odeset('AbsTol',1e-6,'RelTol',1e-3);
tic
[t_ode45,y_ode45]=ode45(@(t,y) fun_f(t,y,par),time_val,y0,options);
t1=toc;

%print CPU time
disp(['time ODE45 ',num2str(t1)])

y_ode45=y_ode45';

%%FE,BE,C-N SOLVER
Ns=0;
N=4;
tol=1.0e-9;
fprintf('FE\n');
[t_fe,y_fe,e_fe,r_fe,t_fe1,y_fe1]=fe_scheme(t0,y0,@(t,y)fun_f(t,y,par),t0,T,Ns,N);
fprintf('-----------\n');
fprintf('BE\n');
[t_be,y_be,e_be,r_be,t_be1,y_be1]=be_scheme(t0,y0,@(t,y)fun_f(t,y,par),t0,T,Ns,N,tol);
fprintf('-----------\n');
fprintf('CN\n');
[t_cn,y_cn,e_cn,r_cn,t_cn1,y_cn1]=cn_scheme(t0,y0,@(t,y)fun_f(t,y,par),t0,T,Ns,N,tol);
fprintf('------------\n');

%% figure showing the state variables
fig1=figure(1);
set(gcf,'color','w');
for jj=1:size(y_ode45,1)
    subplot(3,2,jj)
    plot( time_val,y_ode45(jj,:),'-r' )
    hold on
    plot( (t_fe1'), y_fe1(jj,:),'-b' )
    plot( (t_be1'), y_be1(jj,:),'-g' )
    plot( (t_cn1'), y_cn1(jj,:),'-k' )
    if strcmp(varname(jj),'D')
        plot(t,data_dead,'o')
        legend('ODE45','FE','BE','C-N','Deaths','Location','NorthWest')
    elseif strcmp(varname(jj),'R')
        plot(t,data_recovered,'o')
        legend('ODE45','FE','BE','C-N','Recovered','Location','NorthWest')
    elseif strcmp(varname(jj),'E')
        plot(t,data_positive,'o')
        legend('ODE45','FE','BE','C-N','New Positive','Location','NorthWest')
    else
        legend('ODE45','FE','BE','C-N','Location','NorthWest')
    end
    xlabel('t')
    ylabel(varname{jj})
    set(gca,'XLim',[t0 T])
    datetick('x','dd mmm','keeplimits','keepticks')
end

fig2=figure(2);
set(gcf,'color','w');
for jj=1:size(y_ode45,1)
    subplot(3,2,jj)
    plot( t_cn, log(y_cn(jj,:)),'-k' )
    hold on
    if strcmp(varname(jj),'D')
        plot(t,log(data_dead),'o')
        legend('C-N','Deaths','Location','NorthWest')
    elseif strcmp(varname(jj),'R')
        plot(t,log(data_recovered),'o')
        legend('C-N','Recovered','Location','NorthWest')
    elseif strcmp(varname(jj),'E')
        plot(t,log(data_positive),'o')
        legend('C-N','New Positive','Location','NorthWest')
    else
        legend('C-N','Location','NorthWest')
    end

    xlabel('t')
    ylabel(varname{jj})
    set(gca,'XLim',[t0 T])
    datetick('x','dd mmm','keeplimits','keepticks')
end


figure
subplot(2,1,1)
plot(1:N, log(e_fe(1:4)),'-g')
hold on
plot(1:N, log(e_be(1:4)),'-b')
plot(1:N, log(e_cn(1:4)),'-m')
%plot([1:N], e_ode45,'-r')
legend('e_{fe}','e_{be}','e_{cn}')
title('Error e^{N}')
xlabel('N')
ylabel('e^{N}')

subplot(2,1,2);
plot(1:N, log(r_fe(2:5)),'-g')
hold on
plot(1:N, log(r_be(2:5)),'-b')
plot(1:N, log(r_cn(2:5)),'-m')
%plot([1:N], r_ode45,'-r')
legend('r_{fe}','r_{be}','r_{cn}')
title('Ratio e^{N}/e^{N-1} \tau = 1.0e-9')
xlabel('N')
ylabel('r^{N}')

%saveas(fig1,['plot',num2str(1),'.png'])

%
%---------------------------------------------------
%
function ydot=fun_f(t,y,par)
% function for the SEAIRD model equations
ydot=0*y;
S=y(1);
E=y(2);
A=y(3);
I=y(4);
R=y(5);
D=y(6);

Np=par(1);
beta0=par(2);
delta=par(3);
sigma=par(4);
gammaA=par(5);
gammaI=par(6);
alpha=par(7);
betaA=par(8);

% if t>datenum('02/24/2020')
%     beta0=beta0*0.9;
% elseif t>datenum('03/08/2020')
%     beta0=beta0*0.9*0.66;
% elseif t>datenum('03/22/2020')
%     beta0=beta0*0.9*0.66*0.6;
% elseif t>datenum('05/04/2020')
%     beta0=beta0*0.9*0.66*0.6*2;
% end

foi=beta0*(I+betaA*A)/(Np-D); % force of infection

Sdot=-foi*S;
Edot=foi*S-delta*E;
Adot=(1-sigma)*delta*E-gammaA*A;
Idot=sigma*delta*E-alpha*I-gammaI*I;
Rdot=gammaI*I+gammaA*A;
Ddot=alpha*I;

ydot(1)=Sdot;
ydot(2)=Edot;
ydot(3)=Adot;
ydot(4)=Idot;
ydot(5)=Rdot;
ydot(6)=Ddot;

end
