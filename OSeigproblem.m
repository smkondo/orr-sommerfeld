% OSeigproblem.m
% 
% DESCRIPTION
% Program to compute the Orr-Sommerfeld matrix for three
% dimensional Poiseuille flows using orthogonal Chebyshev polynomials

% Comparing with results found in Schmid and Henningson (2001)
% Eigenvalue spectra: Fig 3.1a pg 64 
% Eigenfunction: Fig. 3.2a,c,e pg 65

% INPUT
% nc = number of chebyshev modes
% nf = number of finite difference modes
% Re = Reynolds number
% kx = streamwise wave number
% kz = spanwise wave number
%
% OUTPUT
% d = 3D Orr-Sommerfeld matrix

%% Part 1: Compute the eigenspectra for kx = 1, kz = 0
% input data
nc = 200; %number of chebyshev points
Re = 10000; %Reynolds Number
kx = 1; %input alpha value
kz = 0; %input for beta value

%% Chebyshev discretization (Reddy)
% generate Chebyshev differentiation matrices
[D0,D1,D2,D4]=Dmat(nc);
% set up Orr-Sommerfeld matrices A and B
[A,B]=pois(nc,kx,kz,Re,D0,D1,D2,D4);

% solve the eigenvalue problem 
[~,echeb] = eigs(A,B,200,'smallestabs');

%% plot
% plot the eigenspectrum
figure(1)
plot(diag(real(echeb)),diag(imag(echeb)),'o','MarkerSize',6)
ylim([-1 0.1]);
xlim([0 1]);
xlabel('Cr', 'FontSize', 18)
ylabel('Ci','FontSize', 18)
title('Eigenspectrum of Orr-Sommerfeld Problem', 'FontSize', 14)
legend('Chebyshev (Reddy)','location','southwest')

%% Part 2: Computing the eigenfunction
% input data
nc = 100; 
Re =  5000; %Reynolds Number
kx = 1; %input alpha value
kz = 1; %input for beta value

%% Chebyshev discretization (Reddy)
% generate Chebyshev differentiation matrices
[D0,D1,D2,D4]=Dmat(nc);
% set up Orr-Sommerfeld matrices A and B
[A,B]=pois(nc,kx,kz,Re,D0,D1,D2,D4);

% solve the eigenvalue problem 
[xcheb,echeb] = eigs(A,B,200,'smallestabs');

vec=(0:nc)';
yj = cos(pi*vec/nc);
%% 
efa1 = cheb_expansion_soln(yj,xcheb(1:101,9));
efp1 = cheb_expansion_soln(yj,xcheb(1:101,69));
efs1 = cheb_expansion_soln(yj,xcheb(1:101,53));

%% plot 
figure(2)
plot(diag(real(echeb)),diag(imag(echeb)),'o','MarkerSize',6)
ylim([-1 0.1]);
xlim([0 1]);
xlabel('Cr', 'FontSize', 18)
ylabel('Ci','FontSize', 18)
title('Eigenspectrum of Orr-Sommerfeld Problem', 'FontSize', 14)
legend('Chebyshev (Reddy)','location','southwest')
%% 
figure(3)
plot(yj,abs(efa1)/max(abs(efa1)),'k-','LineWidth',3)
title('Eigenfunction on A Branch (wall mode)', 'FontSize', 14)
xlabel('y')
ylabel('u')
legend('Chebyshev (Reddy)','location','southwest')

figure(4)
plot(yj,abs(efp1)/max(abs(efp1)),'k-','LineWidth',3)
title('Eigenfunction on P Branch (center mode)', 'FontSize', 14)
xlabel('y')
ylabel('u')
legend('Chebyshev (Reddy)','location','southwest')

figure(5)
plot(yj,abs(efs1)/max(abs(efs1)),'k-','LineWidth',3)
title('Eigenfunction on P Branch (center mode)', 'FontSize', 14)
xlabel('y')
ylabel('u')
legend('Chebyshev (Reddy)','location','southwest')