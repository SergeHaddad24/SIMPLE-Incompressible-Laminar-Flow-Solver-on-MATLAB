function [Fe,Fw,Fn,Fs] = rie_chow(area,u,v,p,ac_u,ac_v,m,n,rho,Sex,Sey,Swx,Swy,Snx,Sny,Ssx,Ssy,mag_ce,mag_ce_up)
Fe = zeros(m-1,n-1);
Fw= zeros(m-1,n-1);
Fs=zeros(m-1,n-1);
Fn=zeros(m-1,n-1);
eps0 =1e-14;
delta_x = 2/(m-1);
delta_y = 1/(n-1);
velocity_entry = round(0.2/delta_x);
pressure_exit = round(0.1/delta_y);
p_out=1e5;
% Clamp d_face = area/ac to prevent explosion at near-stagnation cells.
% Upper bound: area/(4*mu) is the pure-diffusion limit; allow 10x margin.
d_face_max = 10 * max(delta_x,delta_y) / (4 * 5e-5 + eps0);
dpdx=zeros(m-1,n-1);
dpdy=zeros(m-1,n-1);
for i=1:m-1
    for j=1:n-1
        if i==1
            dpdx(i,j) = (p(i+1,j)-p(i,j))/delta_x;
        elseif i==m-1
            dpdx(i,j) = (p(i,j)-p(i-1,j))/delta_x;
        else
            dpdx(i,j)= (p(i+1,j)-p(i-1,j))/(2*delta_x);
        end
        if j==1
            dpdy(i,j)= (p(i,j+1)-p(i,j))/delta_y;
        elseif j==n-1
            dpdy(i,j) = (p(i,j)-p(i,j-1))/delta_y;
        else
            dpdy(i,j) = (p(i,j+1)-p(i,j-1))/(2*delta_y);
        end
    end
end
for i=1:m-1
    for j=1:n-1
        if i<m-1
            uf_bar = 0.5*(u(i,j)+u(i+1,j));
            vf_bar = 0.5*(v(i,j)+v(i+1,j));
            dP = min(area(i,j)/max(abs(ac_u(i,j)),eps0), d_face_max);
            dE = min(area(i+1,j)/max(abs(ac_u(i+1,j)),eps0), d_face_max);
            d_face = 0.5*(dE+dP);
            direct_grad= (p(i+1,j)-p(i,j))/mag_ce(i+1,j);
            interp_grad=0.5*(dpdx(i,j)+dpdx(i+1,j));
            uf= uf_bar -d_face*(direct_grad-interp_grad);
            vf=vf_bar;
            Fe(i,j)=rho*(uf*Sex(i,j)+vf*Sey(i,j));
        else
            Fe(i,j)=0;
        end
        if i>1
            uf_bar = 0.5*(u(i-1,j)+u(i,j));
            vf_bar = 0.5*(v(i-1,j)+v(i,j));
            dP = min(area(i,j)/max(abs(ac_u(i,j)),eps0), d_face_max);
            dW = min(area(i-1,j)/max(abs(ac_u(i-1,j)),eps0), d_face_max);
            d_face = 0.5*(dW+dP);
            direct_grad= (p(i,j)-p(i-1,j))/mag_ce(i,j);
            interp_grad=0.5*(dpdx(i,j)+dpdx(i-1,j));
            uf= uf_bar -d_face*(direct_grad-interp_grad);
            vf=vf_bar;
            Fw(i,j)=rho*(uf*Swx(i,j)+vf*Swy(i,j));
        else
            if j<= pressure_exit
                uf_bar = u(i,j);
                vf_bar = v(i,j);
                d_face = min(area(i,j)/max(abs(ac_u(i,j)),eps0), d_face_max);
                direct_grad = (p(i,j)-p_out)/mag_ce(i,j);
                interp_grad = dpdx(i,j);
                uf = uf_bar -d_face*(direct_grad-interp_grad);
                vf=vf_bar;
                Fw(i,j)=rho*(uf*Swx(i,j)+vf*Swy(i,j));
            else
                Fw(i,j)=0;
            end
        end
    end
end
for i=1:m-1
    for j=1:n-1
        if j<n-1
            uf_bar = 0.5*(u(i,j)+u(i,j+1));
            vf_bar = 0.5*(v(i,j)+v(i,j+1));
            dP = min(area(i,j)/max(abs(ac_v(i,j)),eps0), d_face_max);
            dN = min(area(i,j+1)/max(abs(ac_v(i,j+1)),eps0), d_face_max);
            d_face = 0.5*(dN+dP);
            direct_grad= (p(i,j+1)-p(i,j))/mag_ce_up(i,j+1);
            interp_grad=0.5*(dpdy(i,j)+dpdy(i,j+1));
            uf= uf_bar;
            vf=vf_bar -d_face*(direct_grad-interp_grad);
            Fn(i,j)=rho*(uf*Snx(i,j)+vf*Sny(i,j));
        else 
            if i>=(m-1)-velocity_entry+1
                uf=0;
                vf=-5;
                Fn(i,j)=rho*(uf*Snx(i,j)+vf*Sny(i,j));
            else
                Fn(i,j)=0;
            end
        end
        if j>1
            uf_bar = 0.5*(u(i,j-1)+u(i,j));
            vf_bar = 0.5*(v(i,j-1)+v(i,j));
            dP = min(area(i,j)/max(abs(ac_v(i,j)),eps0), d_face_max);
            dS = min(area(i,j-1)/max(abs(ac_v(i,j-1)),eps0), d_face_max);
            d_face = 0.5*(dS+dP);
            direct_grad= (p(i,j)-p(i,j-1))/mag_ce_up(i,j);
            interp_grad=0.5*(dpdy(i,j)+dpdy(i,j-1));
            uf= uf_bar;
            vf=vf_bar -d_face*(direct_grad-interp_grad);
            Fs(i,j)=rho*(uf*Ssx(i,j)+vf*Ssy(i,j));
        else
            Fs(i,j)=0;
        end
    end
end
end
               
