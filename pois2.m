%pois3.m
function [A,B,G] = pois(n,kx,kz,Re,D0,D1,D2,D4,u,du)
%
% Function to create Orr-Sommerfeld matrices using Chebyshev
% pseudospectral discretization for plane Poiseuille flow
% profile
%
% n = number of modes
% kx = streamwise wavenumber
% kz = spanwise wavenumber
% Re = Reynolds numbers
% DO = zero'th derivative matrix
% D1 = first derivative matrix
% D2 = second derivative matrix
% D3 = third derivative matrix
% D4 = fourth derivative matrix

% mean velocity
ak2=kx^2+kz^2;
N=n+1;
vec=(0:1:n)';
% u=(ones(length(vec),1)-cos(pi*vec/n).^2);
% du=-2*cos(pi*vec/n);

%A11 = Orr-Sommerfeld equation 
%A21 = cross-term
%A22 = Squire equation

%B11 = weight matrix

% set up Orr-Sommerfeld matrix
B11=D2-ak2*D0; 
A11=-(D4-2*ak2*D2+(ak2^2)*D0)/(1i*Re); 
A11=A11+kx*(u*ones(1,length(u))).*B11+kx*2*D0; 

%Applying the BC, which are implemented on the first and last two rows of
%the B matrix. 
%The same rows in matrix A are arbitrary complex multiples of those used in
%matrix B. 
% er=-2000*1i; %arbitrary complex multiple
er = -1000*1i;
B11=[D0(1,:); D1(1,:); B11(3:N-2,:); D1(N,:); D0(N,:)]; %applying the BC on the fist and last two rows
A11=[er*D0(1,:); er*D1(1,:); A11(3:N-2,:); er*D1(N,:); er*D0(N,:) ]; %complex multiples of the same rows 


% set up Squire matrix and cross-term matrix
A21=kz*(du*ones(1,length(u))).*D0(1:N,:); %cross-term
A22=kx*(u*ones(1,length(u))).*D0-(D2-ak2*D0)/(1i*Re); %Squire equation
B22=D0; %1

A22=[er*D0(1,:); A22(2:N-1,:); er*D0(N,:)]; %Applying the same boundary condition on Matrix A
A21=[zeros(1,N); A21(2:N-1,:); zeros(1,N)]; %Zero boundary condition at the wall
% combine all the blocks
A=[A11 zeros(N,N); A21 A22];
B=[B11 zeros(N,N); zeros(N,N) B22];

%set up the C matrix
Der1 = Der(n);
C11 = 1i*kx*Der1;C12 = -1i*kz*eye(size(D0));
C21 = ak2*eye(size(D0)); C22 = zeros(size(D0));
C31 = 1i*kz*Der1;C32 = 1i*kx*eye(size(D0));

C = [C11 C12; C21 C22; C31 C32];

G = 1/ak2*C;

end
