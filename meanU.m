%meanU.m
function [u, du] = meanU(Re_tau, u_tau, n)
%this function inputs the friction Reynolds number and returns the
%streamwise mean velocity profile using an empirical law derived by
%Reynolds and Tiederman in 1967 in 
%"Stability of turbulent channel flow with application to Malkus's theory"

%input
%Re_tau: friction Reynolds number
%u_tau: characteristic wall-velocity of the flow
%n: number of Chebyshev points across the wall-normal direction

%output
%u: streamwise mean velocity
%du: discrete difference of mean velocity between adjacent wall-normal
%points

%define constants of model 
alpha = 25.4;
kappa = 0.426;

vec=(0:n)';
yj = cos(pi*vec/n); yp = Re_tau - yj(1:n/2+1)*Re_tau;

%using the Trapezoid method
u = zeros(size(yp)); %up with Driest's damping function
for j=1:length(yp)
    y = yp(1:j);
    lm = (kappa*y/3.*(6-11*(y/Re_tau)+8*(y/Re_tau).^2-2*(y/Re_tau).^3).*(1.-exp(-1*y./alpha))).^2;
    nu_t = sqrt(1+lm);
    us = 2.*(1-y/Re_tau)./(1+nu_t);
    ut = trapz(y,us);
    u(j) = ut;
end

du = zeros(size(yp));
for i=1:length(yp)
    lm = (kappa*yp(i)/3.*(6-11*(yp(i)/Re_tau)+8*(yp(i)/Re_tau).^2-2*(yp(i)/Re_tau).^3).*(1.-exp(-1*yp(i)./alpha))).^2;
    nu_t = sqrt(1+lm);
    dut = 2.*(1-yp(i)/Re_tau)./(1+nu_t);
    du(i) = dut;
end 

u = [u; flip(u(1:n/2))]*u_tau;
du =[du; flip(du(1:n/2))]*Re_tau*u_tau;
end 