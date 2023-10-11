function dVdx = d_dx_c(V,dx,dim)
%% Description
% d_dx_c is for derivative along one dimension with constant spacing (c stands for constant dx)
% Centered difference is used at the interior, while varying off-centered
% difference is used at the two boundary of the targeted dimension.
% The output variable dVdx is of the same dimension as V.
% The dx is assumed to be single value, and should constant at the targeted
% dimension.
Ndim = length(size(V));

if dim <= Ndim
    dimsize = size(V);
    Nx = dimsize(dim);
    dimOrder = 1:Ndim;
    
    dimOrder2 = dimOrder;
    dimOrder2(dim) = 1;
    dimOrder2(1) = dim;
    
    V2 = permute(V,dimOrder2);
    dimsize2 = size(V2);
    
    dVdx2 = zeros(dimsize2);
    
    dVdx2(2:Nx-1,:,:,:,:,:) = (V2(3:Nx,:,:,:,:,:) - V2(1:Nx-2,:,:,:,:,:))/(2*dx);
    dx1 = dx;
    dx2 = 2*dx;
    dVdx2(1,:,:,:,:,:) = (V2(2,:,:,:,:,:)*(dx2^2) - V2(3,:,:,:,:,:)*(dx1^2) - V2(1,:,:,:,:,:)*(dx2^2 - dx1^2))/(dx1*dx2*(dx2-dx1));
    dx1 = -2*dx;
    dx2 = -dx;
    dVdx2(Nx,:,:,:,:,:) = (V2(Nx-2,:,:,:,:,:)*(dx2^2) - V2(Nx-1,:,:,:,:,:)*(dx1^2) - V2(Nx,:,:,:,:,:)*(dx2^2 - dx1^2))/(dx1*dx2*(dx2-dx1));
    
    % Permute dVdx2 back to dVdx using the same permute order (which is also the inverse order)
    dVdx = permute(dVdx2,dimOrder2);
    
    
else
    'Error dim is wrong!!!'
end


end




