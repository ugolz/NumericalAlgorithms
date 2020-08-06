close all
clear 

% Parameters for the Cauchy problem
t0= datenum('2/1/2020');  % initial time
T=  t0+100;               % final time

%% parametri SIR

par(1)=10e6; %Np
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
h=0.5;
time_val= t0:h:T;

% MATLAB SOLVER
options=odeset('AbsTol',1e-6,'RelTol',1e-3);
tic
[t_ode45,y_ode45]=ode45(@(t,y) fun_f(t,y,par),time_val,y0,options);
t1=toc;

%print CPU time
disp(['time ODE45 ',num2str(t1)])

y_ode45=y_ode45';

%% figure showing the state variables
fig1=figure(1);
set(gcf,'color','w');
for jj=1:size(y_ode45,1)
    subplot(3,2,jj)
    plot( time_val,y_ode45(jj,:),'-r' )
    if jj==size(y_ode45,1)
        legend('ODE45','Location','NorthWest')
    end
    xlabel('t')
    ylabel(varname{jj})
    set(gca,'XLim',[t0 T])
    datetick('x','dd mmm','keeplimits','keepticks')
end
saveas(fig1,['plot',num2str(1),'.png'])

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
