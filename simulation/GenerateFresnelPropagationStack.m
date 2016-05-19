function [HStack]=GenerateFresnelPropagationStack(Nx,Ny,z, lambda, ps)
% Lambda, ps, z has unit meter.
% z is the distance from focal plane.
% lambda is wavelength
% ps is pixel size.

cx=[1:Nx] - (floor(Nx/2)+1);
cy=[1:Ny] - (floor(Ny/2)+1);



[us, vs]=ndgrid(cx, cy);
us=us/Nx/ps; vs=vs/Ny/ps;
%us=ifftshift(us); vs=ifftshift(vs);

HStack = exp(1i*lambda*pi*(us.^2+vs.^2)* z);

%HStack(:,:)= exp(1i*lambda*pi*(us.^2+vs.^2)* z);

end

