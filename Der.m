function Der1=Der(N)
nos = N+1;
vec=(0:N)';
yj = cos(pi*vec/N);
Der1 = zeros(N+1);

Der1(1,1) = (2*N^2+1)/6;
Der1(nos,nos) = -1*(2*N^2+1)/6;
Der1(1,nos) = (-1)^N/(yj(1)-yj(nos));
Der1(nos,1) = (-1)^N/(yj(nos)-yj(1));

for j = 1:N-1
    Der1(1,j+1) = 2*(-1)^j/(yj(1)-yj(j+1)); 
    Der1(nos,j+1) = 2*(-1)^(N+j)/(yj(nos)-yj(j+1));
end

for i = 1:N-1
    Der1(i+1,1) = (-1)^i/(2*(yj(i+1)-yj(1)));
    Der1(i+1,nos) = (-1)^(i+N)/(2*(yj(i+1)-yj(nos)));
end 

for i = 1:N-1
    for j = 1:N-1
        if i == j
            Der1(i+1,j+1) = -1*yj(i+1)/(2*(1-yj(i+1)^2));
        else
            Der1(i+1,j+1) = (-1)^(i+j)/(yj(i+1)-yj(j+1));
        end 
    end 
end 

end 
            
            
        
        

