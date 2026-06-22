function [u_new,v_new,p_new,Fe_new,Fw_new,Fn_new,Fs_new]=correct_fields(u_star,v_star,p_old,p_prime,Fe_star,Fw_star,Fn_star,Fs_star,ac_u,ac_v,rho,area,m,n,alpha_p,Sex,Sey,Swx,Swy,Snx,Sny,Ssx,Ssy,mag_ce,mag_ce_up)
eps0=1e-14;
u_new =u_star;
v_new =v_star;
p_new =p_old+alpha_p*p_prime;
Fe_new=Fe_star;
Fw_new=Fw_star;
Fn_new=Fn_star;
Fs_new=Fs_star;
dx = 2/(m-1);
dy = 1/(n-1);
pressure_exit = round(0.1/dy);
velocity_entry = round(0.2/dx);
% d_face clamp (must match rie_chow clamp)
d_face_max = 10 * max(dx,dy) / (4 * 5e-5 + eps0);

% --- Velocity correction: INTERIOR cells only ---
% Do NOT correct at boundary cells — those are either prescribed (Dirichlet)
% or handled via flux corrections. Correcting them with d_u=area/ac (which
% is huge at low-diffusion boundary cells) is the primary cause of divergence.
for i = 2:m-2
    for j = 2:n-2
        d_u = min(area(i,j)/max(abs(ac_u(i,j)),eps0), d_face_max);
        d_v = min(area(i,j)/max(abs(ac_v(i,j)),eps0), d_face_max);
        gradp_x = (p_prime(i+1,j) - p_prime(i-1,j)) / (2*dx);
        gradp_y = (p_prime(i,j+1) - p_prime(i,j-1)) / (2*dy);
        u_new(i,j) = u_star(i,j) - d_u*gradp_x;
        v_new(i,j) = v_star(i,j) - d_v*gradp_y;
    end
end

% --- Face flux corrections ---
for i=1:m-1
    for j=1:n-1
        Ae=hypot(Sex(i,j),Sey(i,j));
        Aw=hypot(Swx(i,j),Swy(i,j));
        An=hypot(Snx(i,j),Sny(i,j));
        As=hypot(Ssx(i,j),Ssy(i,j));
        if i<m-1
            dP=min(area(i,j)/max(abs(ac_u(i,j)),eps0), d_face_max);
            dE=min(area(i+1,j)/max(abs(ac_u(i+1,j)),eps0), d_face_max);
            d_face=0.5*(dP+dE);
            Ce=rho*d_face*Ae/mag_ce(i+1,j);
            Fe_new(i,j)=Fe_star(i,j)-Ce*(p_prime(i+1,j)-p_prime(i,j));
        end
        if i>1
            dP=min(area(i,j)/max(abs(ac_u(i,j)),eps0), d_face_max);
            dW=min(area(i-1,j)/max(abs(ac_u(i-1,j)),eps0), d_face_max);
            d_face=0.5*(dP+dW);
            Cw=rho*d_face*Aw/mag_ce(i,j);
            Fw_new(i,j)=Fw_star(i,j)+Cw*(p_prime(i,j)-p_prime(i-1,j));
        end
        if i == 1 && j <= pressure_exit
            dP = min(area(i,j)/max(abs(ac_u(i,j)),eps0), d_face_max);
            Cw = rho*dP*Aw/mag_ce(i,j);
            Fw_new(i,j) = Fw_star(i,j) + Cw*(p_prime(i,j) - 0);
        end
        if j<n-1
            dP=min(area(i,j)/max(abs(ac_v(i,j)),eps0), d_face_max);
            dN=min(area(i,j+1)/max(abs(ac_v(i,j+1)),eps0), d_face_max);
            d_face=0.5*(dP+dN);
            Cn=rho*d_face*An/mag_ce_up(i,j+1);
            Fn_new(i,j)=Fn_star(i,j)-Cn*(p_prime(i,j+1)-p_prime(i,j));
        end
        if j>1
            dP=min(area(i,j)/max(abs(ac_v(i,j)),eps0), d_face_max);
            dS=min(area(i,j-1)/max(abs(ac_v(i,j-1)),eps0), d_face_max);
            d_face=0.5*(dP+dS);
            Cs=rho*d_face*As/mag_ce_up(i,j);
            Fs_new(i,j)=Fs_star(i,j)+Cs*(p_prime(i,j)-p_prime(i,j-1));
        end
    end
end

end
