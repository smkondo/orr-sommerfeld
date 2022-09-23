% Function to convert Fourier modes to physical coordinates
% function [U]=fourier2physical(u(r),k,z,n,theta,[omega,t]);
function [eta]=fourier2physical2(u,kx,kz,xval,zval)
% U is the fourier-transformed velocity field
% u is the input (complex) mode shape
% k,n : streamwise, azimuthal wavenumber
% x,theta: streamwise, azimuthal coordinates

% Optional can also include frequency (omega) and time (t) as inputs

n = length(u);

Up = zeros(length(zval),length(xval),n);    % physical value of u

for indx = 1:length(xval)
    x = xval(indx);
    for indz = 1:length(zval)
        z = zval(indz);
        Up(indz,indx,:) =  squeeze(Up(indz,indx,:)) + ...
        u*exp(1i*kx*x + 1i*kz*z);
    end
end
eta = real(Up); % only real part exist

end