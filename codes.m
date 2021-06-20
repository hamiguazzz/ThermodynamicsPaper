clear,close all
tic

rng(0)

n=10^6;
t_end=3;
t_0=0.1;
delta_tao=0.005;
gamma=1;
D=1;
beta=1;
x_cut=20;
tao=linspace(0,t_end,t_end/delta_tao)';
X=0;

F = @(tao,X) -V_differential(X,tao,t_end,t_0)/gamma;
G = @(tao,X) sqrt(2*D);

SDE = sde(F, G, 'StartState', X);
[S,T] = simByEuler(SDE, t_end/delta_tao, 'DeltaTime', delta_tao,'nTrials',n);

figure(1)
for i = 1:5
    plot(T(:),S(:,1,i))
    hold on
end
title('motion of five particles','fontsize',9,'fontname','Times New Roman')
box off
xlabel('t/s','fontsize',9,'fontname','Times New Roman')
ylabel('Position','fontsize',9,'fontname','Times New Roman')
set(gcf,'units','centimeters','position',[10,5,7.59,5.7])
set(gca,'linewidth',0.3)
saveas(gcf,'1.jpg')

figure(2)
[xi,yi]=accurate_distribution(S(t_0/delta_tao,1,:),50);
plot(xi,yi)
hold on
coordinate_x=linspace(min(S(t_0/delta_tao,1,:)),max(S(t_0/delta_tao,1,:)),100);
rho=exp(-coordinate_x.^2/4/D/t_0)/sqrt(4*pi*D*t_0);
plot(coordinate_x,rho)
legend('simulation','theory')
title('diffusion motion','fontsize',9,'fontname','Times New Roman')
legend({'simulation','theory'},'fontsize',9,'fontname','Times New Roman','location','northwest')
legend('boxoff')
box off
xlabel('x','fontsize',9,'fontname','Times New Roman')
ylabel('Density','fontsize',9,'fontname','Times New Roman')
set(gcf,'units','centimeters','position',[10,5,7.59,5.7])
set(gca,'linewidth',0.3)
saveas(gcf,'2.jpg')

figure(3)
[xi,yi]=accurate_distribution(S(end,1,:),50);
plot(xi,yi)
hold on
t=t_end;
coefficient=integral(@(x) distribution(x,D,t_end,beta,t,t_0),-x_cut,x_cut);
coordinate_x=linspace(min(S(t/delta_tao,1,:)),max(S(t/delta_tao,1,:)),100);
distribution_norm=distribution(coordinate_x,D,t_end,beta,t,t_0)/coefficient;
plot(coordinate_x,distribution_norm)
title('diffusion motion in a noncofined potential','fontsize',9,'fontname','Times New Roman')
legend({'simulation','theory'},'fontsize',9,'fontname','Times New Roman','location','northwest')
legend('boxoff')
box off
xlabel('x','fontsize',9,'fontname','Times New Roman')
ylabel('Density','fontsize',9,'fontname','Times New Roman')
set(gcf,'units','centimeters','position',[10,5,7.59,5.7])
set(gca,'linewidth',0.3)
saveas(gcf,'3.jpg')
toc

function V_diff=V_differential(x,tao,t,t_0)
V_diff=A_tao(tao,t,t_0).*exp(-(x-B_tao(tao,t,t_0)).^2/2).*(x-B_tao(tao,t,t_0));
end

function V=gaussian_potential(x,tao,t,t_0)
%定义势垒
V=-A_tao(tao,t,t_0).*exp(-(x-B_tao(tao,t,t_0)).^2/2);
end

function A=A_tao(tao,t,t_0)
A=theta_tao(tao,t_0);
% A=theta_tao(tao,t_0).*sin(pi*(tao-t_0)/(t-t_0));
end

function theta=theta_tao(tao,t_0)
if tao>=t_0
    theta=1;
else
    theta=0;
end

end

function B=B_tao(tao,t,t_0)
B=1;
end

function P_GB=distribution(x,D,t,beta,tao,t_0)
P_GB=exp(-x.^2/4/D/t-beta*gaussian_potential(x,tao,t,t_0));
end

function [xi,yi]=accurate_distribution(v,interval)

accurate=5;
v=sort(v);
len=(v(end)-v(1))/interval;

n=length(v);
xi=zeros((interval-1)*accurate+1,1);
yi=zeros((interval-1)*accurate+1,1);
i=1;
j=1;
for k=1:(interval-1)*accurate+1
    front=v(1)+(k-1)*len/accurate;
    rear=v(1)+len*(1+(k-1)/accurate);
    xi(k)=(front+rear)/2;
    while j<n && v(j+1)<=rear
        j=j+1;
    end
    while i<j && v(i)<front
        i=i+1;
    end
    yi(k)=j-i+1;
end

yi=yi./n/len;

end
