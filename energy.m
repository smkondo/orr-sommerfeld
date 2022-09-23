%energy.m
%calculating the weights for the energy norm
function M = energy(n,ak2)
    Cos=two(n+1); 
    Dos=deven(n+1); 
    Wos=(Dos'*Cos*Dos).*(1/ak2)+Cos;
    Wsq=1/ak2.*Cos;

    F = [Wos zeros(n+1); zeros(n+1) Wsq];
    M = chol(F);
end 