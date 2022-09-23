% test.m
% Computing and plotting the resolvent modes for Channel flow Re_tau = 182
%1. The mean velocity profile from DNS data and empirical formula
%2. The 20 singular values from resolvent analysis (DNS)
%3. The resolvent modes (all three values on one plot)
%4. The resolvent modes contour plot

%% 1.Compute the velocity profile of Re_tau = 182
%define the parameters
n=100;
nu = 3.50000e-04; 
u_tau = 6.37309e-02;
Re_tau = u_tau*1/nu;
Re = 1/nu;

vec=(0:n)'; yj = cos(pi*vec/n); yp = Re_tau - yj(1:n/2+1)*Re_tau;

%mean velocity from empirical formula
[u, du] = meanU(Re_tau, u_tau, n);

%mean velocity from DNS data
[udns, dudns] = meanUDNS(Re_tau,u_tau,n,yj);

%% 1.1 plot the mean velocity profile
figure(1)
yplot = abs(yj*Re_tau-Re_tau);
semilogx(yplot,u,'LineWidth',2)
hold on
semilogy(yplot,udns,'LineWidth',2)
hold off
title('mean velocity profile Re_\tau = 180')
ylabel('u^+')
xlabel('y^+')
legend('empirical','DNS (Lee and Moser)','location','northwest')
xlim([0 180])

%% 2.compute the 20 largest singular values
kx = 2*pi/1000*Re_tau; %input alpha value (lambda_x^+ = 1000)
kz = 2*pi/100*Re_tau; %input for beta value (lambda_z^+ = 100)
c = 10; om = kx*c*u_tau; %c = phase speed, om = frequency

ak2=kx^2+kz^2;

[D0,D1,D2,D4]=Dmat(n);

Cos=two(n+1); 
Dos=deven(n+1); 
Wos=Dos'*Cos*Dos+ak2*Cos;
Wsq=two(n+1);

F = [Wos zeros(n+1); zeros(n+1) Wsq];
M = chol(F);

%SVD of resolvent operator (DNS)
[A,B,C]=pois(n,kx,kz,Re,D0,D1,D2,D4,udns,dudns);

RA = M/(om*eye(2*n+2)-B\A)/M;
[~,ss,~] = svds(RA,20);

%% 2.1 plot the 20 largest singular values
figure(2)
semilogy(diag(ss),'o')
title('20 largest \sigma_i from SVD of resovelnt operator (DNS)')
xlabel('i')
ylabel('\sigma_i')


%% 3.compute the resolvent quantity of the fluctuating velocity
velfluc = resolvent(n,kx,kz,ak2,om,Re,D0,D1,D2,D4,udns,dudns,yj);

normalu = velfluc(1:n+1);
normalv = velfluc(n+2:2*n+2); 
normalw = velfluc(2*n+3:3*n+3);

%% 3.1 plot the fourier coefficient of the fluctuating velocity
yplot = abs(yj*Re_tau-Re_tau);
figure(3)
plot(yplot,abs(normalu),'LineWidth',2)
hold on
plot(yplot,abs(normalv),'LineWidth',2)
plot(yplot,abs(normalw),'LineWidth',2)
hold off
title('first resolvent mode of velocity fluctuations')
xlabel('y^+')
ylabel('response')
xlim([0 180])
legend('u''', 'v''', 'w''')


%% 4.Contour plot of fluctuating quantity in physical space
%Fourier space to physical space
xval = linspace(0,1000/Re_tau,100); 
zval = linspace(0,1,100); %Nz = 520 nodes

[etau] = fourier2physical(normalu,kx,kz,xval,zval); 
[etav] = fourier2physical(normalv,kx,kz,xval,zval);
[etaw] = fourier2physical(normalw,kx,kz,xval,zval);

%calculating the derivatives of u,v,w with respect to x,y,z
etau2 = squeeze(etau(:,50,:));

dudz = 1/100.*(etau2(2:100,:) - etau2(1:99,:)); %done
dudz2 = [1/100.*etau2(1,:); dudz];
dwdx = 1000/Re_tau/100.*(etaw(:,2:100,:) - etaw(:,1:99,:));

dvdx = 1000/Re_tau/100.*(etav(:,2:100,:) - etav(:,1:99,:));
dudy = zeros(size(etau2)); %done
for i = 2:n+1
    dudy(:,i-1)=(etau2(:,i)-etau2(:,i))/abs(yj(i)-yj(i-1));
end

omegav = dudz2 - squeeze(dwdx(:,50,:));
omegaw = squeeze(dvdx(:,50,:)) - dudy;
%Use Fourier coefficients and transform to physical space 
xval = 1500/Re_tau; %at a specific x location
zval = linspace(0,1,100);

[etau] = fourier2physical(normalu,kx,kz,xval,zval); 
etau = squeeze(etau); etau = flip(permute(etau,[2,1]),2);

[etav] = fourier2physical(normalv,kx,kz,xval,zval);
etav = squeeze(etav); etav = flip(permute(etav,[2,1]),2);

[etaw] = fourier2physical(normalw,kx,kz,xval,zval);
etaw = squeeze(etaw); etaw = flip(permute(etaw,[2,1]),2);

%% 4.1 plot the streamwise velocity fluctuation contour
zplot = zval*Re_tau; yplot = abs(yj*Re_tau-Re_tau); 
etauplot = etau2';
omv = omegav';
omw = omegaw';
[A,B] = meshgrid(zplot,yplot);
figure(4)
contourf(A,B,etauplot)
hold on
quiver(A(1:5:end),B(1:5:end),omw(1:5:end),omv(1:5:end),1.5,'k','LineWidth',1)
hold off
title('Contour plot of u''')
xlim([0 150])
ylim([0 150])
xlabel('z^+')
ylabel('y^+')
colorbar
colormap(redblue(50))
