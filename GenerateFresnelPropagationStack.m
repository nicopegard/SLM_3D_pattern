function [HStack]=GenerateFresnelPropagationStack(Nx,Ny,z, lambda, ps, f)
% Lambda, ps, z has unit meter.
% z is the distance from focal plane.
% lambda is wavelength
% ps is pixel size.

cx=[1:Nx] - (floor(Nx/2)+1);
cy=[1:Ny] - (floor(Ny/2)+1);



[us, vs]=ndgrid(cx, cy);
us = (us * ps)./sqrt((us * ps).^2 + f^2);
vs = (vs * ps)./sqrt((vs * ps).^2 + f^2);
us=us / lambda;
vs=vs / lambda;
%us=ifftshift(us); vs=ifftshift(vs);

HStack = exp(1i*lambda*pi*(us.^2+vs.^2)* z);

%HStack(:,:)= exp(1i*lambda*pi*(us.^2+vs.^2)* z);

end

