function [du_dx,du_dy,dv_dx,dv_dy] = velocity_gradients(u,v,m,n,dx,dy)

Nx = m-1;
Ny = n-1;

du_dx = zeros(Nx,Ny);
du_dy = zeros(Nx,Ny);
dv_dx = zeros(Nx,Ny);
dv_dy = zeros(Nx,Ny);

for i = 1:Nx
    for j = 1:Ny

        if i == 1
            du_dx(i,j) = (u(i+1,j)-u(i,j))/dx;
            dv_dx(i,j) = (v(i+1,j)-v(i,j))/dx;
        elseif i == Nx
            du_dx(i,j) = (u(i,j)-u(i-1,j))/dx;
            dv_dx(i,j) = (v(i,j)-v(i-1,j))/dx;
        else
            du_dx(i,j) = (u(i+1,j)-u(i-1,j))/(2*dx);
            dv_dx(i,j) = (v(i+1,j)-v(i-1,j))/(2*dx);
        end

        if j == 1
            du_dy(i,j) = (u(i,j+1)-u(i,j))/dy;
            dv_dy(i,j) = (v(i,j+1)-v(i,j))/dy;
        elseif j == Ny
            du_dy(i,j) = (u(i,j)-u(i,j-1))/dy;
            dv_dy(i,j) = (v(i,j)-v(i,j-1))/dy;
        else
            du_dy(i,j) = (u(i,j+1)-u(i,j-1))/(2*dy);
            dv_dy(i,j) = (v(i,j+1)-v(i,j-1))/(2*dy);
        end

    end
end

end