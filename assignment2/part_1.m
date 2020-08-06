close all
clear all

t0=0;
tf=15;
y0=0.1;%initial condition of Cauchy problem
v0=0;
Ns=1;
N=8;
tol=1.0e-9;
k=1;
m=0.25;

%%set condition for Cauchy problem
t_0(1)=t0;
y_0(1,1)=y0;
y_0(2,1)=v0;
%%analytical solution
tall=linspace(t0,tf,(2^N));
ysol=fsol(tall,y0);
fprintf('FE\n');
[t_fe,y_fe,e_fe,r_fe,t_fe1,y_fe1]=fe_scheme(t_0,y_0,@(x,y)spring_eq(x,y),t0,tf,Ns,N);
fprintf("--------------\n");
fprintf('BE\n');
[t_be,y_be,e_be,r_be,t_be1,y_be1]=be_scheme(t_0,y_0,@(x,y)spring_eq(x,y),t0,tf,Ns,N,tol);
fprintf('--------------\n');
fprintf('C-N\n');
[t_cn,y_cn,e_cn,r_cn,t_cn1,y_cn1]=cn_scheme(t_0,y_0,@(x,y)spring_eq(x,y),t0,tf,Ns,N,tol);
fprintf('--------------\n');

%%ode45
t_fe(1)=t0;
y_ode450(1,1)=y0;
y_ode450(1,2)=v0;
for n=1:N
    h=2^(-n);
    options=odeset('AbsTol',1e-8,'RelTol',1e-6);
    [t_ode45,y_ode45]=ode45(@(t,y) spring_eq(t,y),t0:h:tf,y_ode450,options);
    y_err = y_ode45(:,1);
    e_ode45(n)= max(norm(ysol-y_err));
    if n>=2
        r_ode45(n) = e_ode45(n)/e_ode45(n-1);
    end
    if n==1
        t_ode451=t_ode45;
        y_ode451=y_ode45(:,1);
    end
end


figure
subplot(2,1,1);
plot(t_ode451,y_ode451,'-or')
hold on
plot(t_fe1,y_fe1(1,:),'-*g')
plot(t_be1,y_be1(1,:),'-*b')
plot(t_cn1,y_cn1(1,:),'-*m')
plot(tall,ysol,'-k')
legend('ode','fe','be','cn','sol')
title('N=1')
xlabel('t')
ylabel('y(t)')
ylim([-1 1])

subplot(2,1,2);
plot(t_ode45,y_ode45(:,1),'-or')
hold on
plot((t_fe'),y_fe(1,:),'-*g')
plot((t_be'),y_be(1,:),'-*b')
plot((t_cn'),y_cn(1,:),'-*m')
plot(tall,ysol,'-k')
legend('ode','fe','be','cn','sol')
title('N=8')
xlabel('t')
ylabel('y(t)')
%xlim([0 0.07])
%ylim([0.099 0.1002])

figure
subplot(2,1,1);
plot(1:N, log(e_fe(2:9)),'-g')
hold on
plot(1:N, log(e_be(2:9)),'-b')
plot(1:N, log(e_cn(2:9)),'-m')
%plot([1:N], e_ode45,'-r')
legend('e_{fe}','e_{be}','e_{cn}')
title('Error e^{N} \tau = 1.0e-3')
xlabel('N')
ylabel('e^{N}')

subplot(2,1,2);
plot(1:N, log(r_fe(2:9)),'-g')
hold on
plot(1:N, log(r_be(2:9)),'-b')
plot(1:N, log(r_cn(2:9)),'-m')
%plot([1:N], r_ode45,'-r')
legend('r_{fe}','r_{be}','r_{cn}')
title('Ratio e^{N}/e^{N-1} \tau = 1.0e-3')
xlabel('N')
ylabel('r^{N}')
