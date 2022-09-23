%meanUDNS.m
function [udns, dudns] = meanUDNS(Re_tau,u_tau,n,yj)
yp = Re_tau - yj(1:n/2+1)*Re_tau;
%mean velocity from DNS data
data = importdata('0180_mean_prof.txt');
ydns = data(:,2); uDNS = data(:,3)*u_tau; duDNS = data(:,4)*Re_tau*u_tau;

uinterp = interp1(ydns,uDNS,yp); duinterp = interp1(ydns,duDNS,yp);
udns = [0; uinterp(2:n/2); max(uinterp); flip(uinterp(2:n/2)); 0];
dudns = [max(duinterp); duinterp(2:n/2); 0; flip(duinterp(2:n/2)); max(duinterp)];
end