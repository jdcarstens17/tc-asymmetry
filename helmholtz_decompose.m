function [u_psi, v_psi, u_chi, v_chi, psi,chi] = helmholtz_decompose(U,V, dx, dy)

% This subroutine compute the helmholtz decomposition of a 2D wind field (U,V)
% which are assumed to be on the same 2D grid, with grid spacing dx and dy on
% the two direction.
% This subroutine returns the non-divergent and irrotational wind and the
% corresponding streamfunction psi and velocity potential chi.
% The "environmental" component can be compute easily by subtracting
% (u,v)_psi and (u,v)_chi from the full (u,v).

dimsize = size(U);
if dimsize == size(V)
    
    % Compute Divergence and Vorticity at the interior of the domain
    Div = d_dx_c(U, dx, 1) + d_dx_c(V, dy, 2);
    Vor = d_dx_c(V, dx, 1) - d_dx_c(U, dy, 2);
    
    % Define solving variables for LU decomposition
    Nx = dimsize(1);
    Ny = dimsize(2);
    
    f1= zeros(Nx*Ny,1);           % vorticity 
    y1= zeros(Nx*Ny,1);
    s1= zeros(Nx*Ny,1);
    
    f2= zeros(Nx*Ny,1);           % divergence
    y2= zeros(Nx*Ny,1);
    s2= zeros(Nx*Ny,1);
    
    dum_1d = zeros(Nx,1);
    psi= zeros(Nx,Ny);           % Streamfunction 
    chi= zeros(Nx,Ny);           % Velocity Potential 

    u_psi= zeros(Nx,Ny);           % nondivergent wind 
    v_psi= zeros(Nx,Ny);           
    
    u_chi= zeros(Nx,Ny);           % irrotational wind 
    v_chi= zeros(Nx,Ny);          
    
    A= zeros(Nx,Nx,Ny-1);
    B= zeros(Nx,Nx,Ny);
    C= zeros(Nx,Nx,Ny-1);
    
    V= zeros(Nx,Nx,Ny);
    InvV = zeros(Nx,Nx,Ny);
    U= zeros(Nx,Nx,Ny-1);
    
    for j = 1:Ny
        for i = 1:Nx
            ip = i+1;
            in = i-1;
            jp = j+1;
            jn = j-1;
            
            if (i+1 > Nx)
                ip = i;
            end
            if (i-1 < 1)
                in = i;
            end
            
            B(i,ip,j) = 1/dx^2;
            B(i,in,j) = 1/dx^2;
            B(i,i,j) = -2/dx^2 - 2/dy^2;
            
            if (jn >= 1)
                A(i,i,jn) = 1/dy^2;
            end
            if (jp <= Ny)
                C(i,i,jp-1) = 1/dy^2;
            end
            ii = (j-1)*Nx+i;
            f1(ii) = Vor(i,j);
            f2(ii) = Div(i,j);
        end
    end    
    %% LU decomposition
    % First step is the same for both f1 and f2
    V(:,:,1) = B(:,:,1);
    InvV(:,:,1) = inv(V(:,:,1));
    for j = 1:Ny-1
        j;
        U(:,:,j) = InvV(:,:,j)*C(:,:,j);
        V(:,:,j+1) = B(:,:,j+1) - A(:,:,j)*U(:,:,j);
        InvV(:,:,j+1) = inv(V(:,:,j+1));
    end
    
    % Do this two times for f1 and f2 
    % for f1
    y1(1:Nx,1) = InvV(:,:,1)*f1(1:Nx,1);
    y2(1:Nx,1) = InvV(:,:,1)*f2(1:Nx,1);
    for j = 2:Ny
        dum_1d(1:Nx,1) = y1((j-2)*Nx+1:(j-1)*Nx,1);
        y1((j-1)*Nx+1:j*Nx,1) = InvV(:,:,j)*(  f1((j-1)*Nx+1:j*Nx,1) - A(:,:,j-1)*dum_1d );
        
        dum_1d(1:Nx,1) = y2((j-2)*Nx+1:(j-1)*Nx,1);
        y2((j-1)*Nx+1:j*Nx,1) = InvV(:,:,j)*(  f2((j-1)*Nx+1:j*Nx,1) - A(:,:,j-1)*dum_1d );
    end
    
    s1((Ny-1)*Nx+1:Nx*Ny,1) = y1((Ny-1)*Nx+1:Nx*Ny,1);
    s2((Ny-1)*Nx+1:Nx*Ny,1) = y2((Ny-1)*Nx+1:Nx*Ny,1);
    for j = Ny-1:-1:1
        dum_1d(1:Nx,1) = s1(j*Nx+1:(j+1)*Nx,1);
        s1((j-1)*Nx+1:j*Nx,1) = y1((j-1)*Nx+1:j*Nx,1) - U(:,:,j)*dum_1d;
        
        dum_1d(1:Nx,1) = s2(j*Nx+1:(j+1)*Nx,1);
        s2((j-1)*Nx+1:j*Nx,1) = y2((j-1)*Nx+1:j*Nx,1) - U(:,:,j)*dum_1d;
    end
    
    
    for j = 1:Ny
        for i = 1:Nx
            psi(i,j) = s1((j-1)*Nx+i,1);
            chi(i,j) = s2((j-1)*Nx+i,1);
        end
    end
    
    u_chi = d_dx_c(chi, dx, 1);
    v_chi = d_dx_c(chi, dy, 2);
    
    u_psi = -d_dx_c(psi, dy, 2);
    v_psi = d_dx_c(psi, dx, 1);
    
    
else
    'Error: U and V are not the same size!!!'
    
end
