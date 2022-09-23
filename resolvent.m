%resolvent.m
function velfluc = resolvent(n,kx,kz,ak2,om,Re,D0,D1,D2,D4,u,du,yj)
%this function inputs the mean-velocity profile, the wavenumber-frequency
%pair and returns the corresponding Fourier Coefficients of the most
%energetic (first-rank) fluctuating velocities

%input 
%Re: Reynolds number of the flow
%n: number of chebyshev points along the wall-normal location
%kx: streawise wavenumber
%kz: spanwise wavenumber
%ak2: kx^2 + kz^2
%om: streamwise frequency
%D0, D1, D2, D4: Discrete differential operator in the wall-normal
%u: streamwise mean velocity profile
%du: discrete difference of mean velocity in wall-normal
%yj: wall-normal location of the chebyshev points

%energy norm weights
Cos=two(n+1); 
Dos=deven(n+1); 
Wos=Dos'*Cos*Dos+ak2*Cos;
Wsq=two(n+1);

F = [Wos zeros(n+1); zeros(n+1) Wsq];
M = chol(F);

%mass matrix (B)
%linear operator containing OS,SQ and crossterms (A)
%transformation matrix from eta to u' and v' (C)
[A,B,C]=pois2(n,kx,kz,Re,D0,D1,D2,D4,u,du);

%Resolvent operator (RA)
RA = M/(om*eye(2*n+2)-B\A)/M;

%first-rank SVD of the resolvent operator 
[su,ss,~] = svds(RA,1);

%Chebyshev expansions of the first resolvent mode
PrinRes1 = ss*M\su;
normalv = cheb_expansion_soln(yj,PrinRes1(1:n+1));
normaleta = cheb_expansion_soln(yj,PrinRes1(n+2:2*n+1));

%transform vorticity into velocity components
phi = [normalv; normaleta];
velfluc = C*phi; 
end 