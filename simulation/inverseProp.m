function [ phase ] = inverseProp( imageField, Nx, Ny, dz, lambda, ps )
%INVERSEPROP Summary of this function goes here
%   Detailed explanation goes here
    regularizer = 1e-6;
    Hstack = GenerateFresnelPropagationStack(Nx,Ny,dz, lambda, ps);
    objectField = ifft2(ifftshift(imageField))./(Hstack + regularizer);
    phase = angle(objectField);
end

