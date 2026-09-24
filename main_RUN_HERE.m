Nx = input('Enter number of cells in x-direction: ');
Ny = input('Enter number of cells in y-direction: ');
urf = input('Enter under-relaxation factor (e.g. 0.8): ');
tol = input('Enter convergence tolerance (e.g. 1e-6): ');
maxIter = input('Enter maximum number of SIMPLE iterations: ');
m = Nx + 1;
n = Ny + 1;
rho= 0.8;
mu = 5*10^-5;
cp= 1.03;
k_thermal = 0.036;
dx = 2/Nx;
dy = 1/Ny;
pressure_exit=round(0.1/dy);
velocity_entry=round(0.2/dx);
[X,Y]=TFI(m,n);
centroid_X= cell_centroid_X(X,m,n);
centroid_Y= cell_centroid_Y(Y,m,n);
[area] = get_area(X,Y,m,n);
[Sex,Sey, Swx,Swy,Snx,Sny,Ssx,Ssy,nnx,nny,nex,ney] = surface_vectors(X,Y,m,n);
[f_centroid_east_x,f_centroid_east_y,f_centroid_west_x,f_centroid_west_y,f_centroid_top_x,f_centroid_top_y,f_centroid_bottom_x,f_centroid_bottom_y] = face_centroid(X,Y,m,n);
plotgrid(X,Y,centroid_X,centroid_Y,f_centroid_east_x,f_centroid_east_y,f_centroid_west_x,f_centroid_west_y,f_centroid_top_x,f_centroid_top_y,f_centroid_bottom_x,f_centroid_bottom_y);
 [T_X,T_Y,T_up_X,T_up_Y,E_X,E_Y,E_up_X,E_up_Y,e_X,e_Y,e_up_X,e_up_Y,ce_distance_X,mag_ce,mag_ce_up,ce_distance_Y,ce_distance_up_X,ce_distance_up_Y] = ce_distance(m,n,centroid_X,centroid_Y,f_centroid_east_x,f_centroid_east_y,f_centroid_west_x,f_centroid_west_y,f_centroid_top_x,f_centroid_top_y,f_centroid_bottom_x,f_centroid_bottom_y,Swx,Swy,Sex,Sey,Snx,Sny,Ssx,Ssy);

[ge,gw,gs,gn] =get_interpolationFactors(m,n,area);
maxGS = 1000;
maxSIMPLE = maxIter;   % SIMPLE outer iterations (user input)
p_out = 1e5;           % reference outlet pressure

[u_0, v_0, p_0, T_0 ] = initial_conditions(m,n);
u=u_0;
v=v_0;
p=p_0;
T=T_0;
[Fe, Fw, Fn, Fs]=convective_fluxes(rho,m,n,u,v,Sex,Sey,Swx,Swy,Snx,Sny,Ssx,Ssy);
maxOuter=input('max outer iterations:');
for iterSIMPLE = 1:maxSIMPLE
    u_old = u;
    v_old = v;
    p_old = p;
    T_old = T;
    [du_dx,du_dy,dv_dx,dv_dy] = velocity_gradients(u,v,m,n,dx,dy);
    [ae_u,aw_u,as_u,an_u,ac_u,bc_base_u] = coefficients_u(du_dx,dv_dx,mu,p,Sex,Sey,Ssx,Ssy,Snx,Sny,Swx,Swy,E_X,E_Y,mag_ce,mag_ce_up,E_up_X,E_up_Y,area,m,n,Fe,Fw,Fn,Fs);
    [ae_v,aw_v,as_v,an_v,ac_v,bc_base_v] = coefficients_v(du_dy,dv_dy,mu,p,Sex,Sey,Ssx,Ssy,Snx,Sny,Swx,Swy,E_X,E_Y,mag_ce,mag_ce_up,E_up_X,E_up_Y,area,m,n,Fe,Fw,Fn,Fs);
    ac_u_star = ac_u / urf;
    ac_v_star = ac_v / urf;
    u_star=SMART(u_old, ae_u, aw_u, as_u, an_u, ac_u_star, bc_base_u, Fe, Fw, Fn, Fs, m, n, urf, maxGS, maxOuter, tol);
    v_star=SMART(v_old, ae_v, aw_v, as_v, an_v, ac_v_star, bc_base_v, Fe, Fw, Fn, Fs, m, n, urf, maxGS, maxOuter, tol);
    [maxv, idx] = max(abs(v_star(:)));
    [iMax,jMax] = ind2sub(size(v_star),idx);

    fprintf('max v_star = %.6e at i=%d, j=%d\n', maxv, iMax, jMax);
    fprintf('v_star there = %.6e\n', v_star(iMax,jMax));
    fprintf('ac_v there = %.6e\n', ac_v(iMax,jMax));
    fprintf('bc_v there = %.6e\n', bc_base_v(iMax,jMax));
    fprintf('After momentum: max u_star = %.3e, max v_star = %.3e\n', ...
    max(abs(u_star(:))), max(abs(v_star(:))));
    [Fe_star,Fw_star,Fn_star,Fs_star]= rie_chow(area,u_star,v_star,p_old,ac_u_star,ac_v_star,m,n,rho,Sex,Sey,Swx,Swy,Snx,Sny,Ssx,Ssy,mag_ce,mag_ce_up);
    fprintf('After RC: max Fe = %.3e, max Fw = %.3e, max Fn = %.3e, max Fs = %.3e\n', ...
    max(abs(Fe_star(:))), max(abs(Fw_star(:))), ...
    max(abs(Fn_star(:))), max(abs(Fs_star(:))));
    fprintf('min ac_v = %.6e, max ac_v = %.6e\n', min(ac_v(:)), max(ac_v(:)));
    fprintf('max bc_v = %.6e\n', max(abs(bc_base_v(:))));
    Fe= Fe_star;
    Fw=Fw_star;
    Fn=Fn_star;
    Fs=Fs_star;
   [ae_pprime,aw_pprime,as_pprime,an_pprime,ac_pprime,bc_pprime] = coefficients_pprime( ...
    Fe,Fw,Fn,Fs,ac_u_star,ac_v_star,rho,area,m,n, ...
    Sex,Sey,Swx,Swy,Snx,Sny,Ssx,Ssy,mag_ce,mag_ce_up);
    p_prime0=zeros(m-1,n-1);
    p_prime = solve_pressure_correction(p_prime0,ae_pprime,aw_pprime,as_pprime,an_pprime,ac_pprime,bc_pprime,m,n,maxGS,tol);
   
    alpha_p = 1 - urf;
    [u_new,v_new,p_new,Fe_new,Fw_new,Fn_new,Fs_new]=correct_fields(u_star,v_star,p_old,p_prime,Fe_star,Fw_star,Fn_star,Fs_star,ac_u_star,ac_v_star,rho,area,m,n,alpha_p,Sex,Sey,Swx,Swy,Snx,Sny,Ssx,Ssy,mag_ce,mag_ce_up);
    u=u_new;
    v=v_new;
    p=p_new;
    
    p(1, 1:pressure_exit) = p_out;
    Fe=Fe_new;
    Fw=Fw_new;
    Fn=Fn_new;
    Fs=Fs_new;
    Fe(Nx,:) = 0;
    Fs(:,1)  = 0;
    Fw(1,pressure_exit+1:Ny) = 0;

    Fn(:,Ny) = 0;
    for ii = Nx-velocity_entry+1:Nx
        Fn(ii,Ny) = rho*(0*Snx(ii,Ny) + (-5)*Sny(ii,Ny));
    end
    mass_imbalance = Fe + Fw + Fn + Fs;

    cont_res = sum(abs(mass_imbalance(:))) / ...
        (sum(abs(Fe(:))) + sum(abs(Fw(:))) + ...
        sum(abs(Fn(:))) + sum(abs(Fs(:))) + 1e-30);

    cont_res_hist(iterSIMPLE) = cont_res;

    fprintf('SIMPLE iter %d, continuity residual = %.6e\n', ...
        iterSIMPLE, cont_res);

    if cont_res < tol
        break
    end
   mass_res = max(abs(mass_imbalance(:)));
    fprintf('SIMPLE iter %d, mass residual= %.6e\n',iterSIMPLE,mass_res);
    if mass_res<tol
        break
    end
    if max(abs([u_star(:); v_star(:); Fe_star(:); Fw_star(:); Fn_star(:); Fs_star(:)])) > 1e6
        error('Divergence before pressure correction');
    end
end
[ae_T,aw_T,as_T,an_T,ac_T,bc_T] = coefficients_T( ...
    k_thermal,cp,Fe,Fw,Fn,Fs, ...
    Sex,Sey,Swx,Swy,Snx,Sny,Ssx,Ssy, ...
    mag_ce,mag_ce_up,area,m,n);

urf_T = 0.3;
maxGS_T = 10000;
tol_T = 1e-6;

T = solve_GS_scalar(T,ae_T,aw_T,as_T,an_T,ac_T,bc_T, ...
    m,n,urf_T,maxGS_T,tol_T);



figure;
quiver(centroid_X, centroid_Y, u, v);
xlabel('x');
ylabel('y');
title('Velocity vectors');
axis equal tight;
speed = sqrt(u.^2 + v.^2);

figure;
contourf(centroid_X, centroid_Y, speed, 30, 'LineColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title('Velocity magnitude contour');
axis equal tight;
figure;
semilogy(1:iterSIMPLE, cont_res_hist(1:iterSIMPLE), 'o-');
xlabel('SIMPLE iteration');
ylabel('Continuity residual');
title('Continuity Residual Convergence');
grid on;
figure;
contourf(centroid_X,centroid_Y,T,30,'LineColor','none');
colorbar;
xlabel('x');
ylabel('y');
title('Temperature contour');
axis equal tight;
