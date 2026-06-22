function [u,v] = apply_velocity_BC(u,v,m,n)

Nx = m-1;
Ny = n-1;

dx = 2/Nx;
dy=1/Ny;
velocity_entry = round(0.2/dx);
pressure_exit=round(0.1/dy);

% left boundary
i = 1;
for j = pressure_exit+1:Ny
    
    u(i,j) = 0;
    v(i,j) = 0;
end

% right boundary
i = Nx;
for j = 1:Ny
    u(i,j) = 0;
    v(i,j) = 0;
end

% bottom boundary
j = 1;
for i = 1:Nx
    u(i,j) = 0;
    v(i,j) = 0;
end

% top boundary
j = Ny;
for i = 1:Nx
    if i >= Nx-velocity_entry+1
        u(i,j) = 0;
        v(i,j) = -5;
    else
        u(i,j) = 0;
        v(i,j) = 0;
    end
end

end