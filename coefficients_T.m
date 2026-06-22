function [ae,aw,as,an,ac,bc] = coefficients_T( ...
    k,cp,Fe,Fw,Fn,Fs, ...
    Sex,Sey,Swx,Swy,Snx,Sny,Ssx,Ssy, ...
    mag_ce,mag_ce_up,area,m,n)

Nx = m-1;
Ny = n-1;

ae = zeros(Nx,Ny);
aw = zeros(Nx,Ny);
as = zeros(Nx,Ny);
an = zeros(Nx,Ny);
ac = zeros(Nx,Ny);
bc = zeros(Nx,Ny);

dx = 2/Nx;
dy = 1/Ny;

velocity_entry = round(0.2/dx);
pressure_exit  = round(0.1/dy);

T_bottom = 500;   % bottom wall
T_top    = 300;   % top wall
T_inlet  = 300;   % top-right inlet

h_inf = 15;       % west convective wall
T_inf = 300;

for i = 1:Nx
    for j = 1:Ny

        ac_extra = 0;

        Ae = hypot(Sex(i,j),Sey(i,j));
        Aw = hypot(Swx(i,j),Swy(i,j));
        An = hypot(Snx(i,j),Sny(i,j));
        As = hypot(Ssx(i,j),Ssy(i,j));

        De = k*Ae/mag_ce(i+1,j);
        Dw = k*Aw/mag_ce(i,j);
        Dn = k*An/mag_ce_up(i,j+1);
        Ds = k*As/mag_ce_up(i,j);

        FeT = cp*Fe(i,j);
        FwT = cp*Fw(i,j);
        FnT = cp*Fn(i,j);
        FsT = cp*Fs(i,j);

        ae(i,j) = De + max(-FeT,0);
        aw(i,j) = Dw + max(-FwT,0);
        an(i,j) = Dn + max(-FnT,0);
        as(i,j) = Ds + max(-FsT,0);

        % -----------------------------
        % West boundary
        % lower-left outlet: zero-gradient
        % upper-left wall: convection to ambient
        % -----------------------------
        if i == 1
            if j <= pressure_exit
                % outlet, dT/dn = 0
                aw(i,j) = 0;
            else
                % convection: -k dT/dn = h(T - Tinf)
                aB = h_inf*Aw;

                bc(i,j) = bc(i,j) + aB*T_inf;
                ac_extra = ac_extra + aB;

                aw(i,j) = 0;
            end
        end

        % -----------------------------
        % East/right wall: insulated
        % -----------------------------
        if i == Nx
            ae(i,j) = 0;
        end

        % -----------------------------
        % South/bottom wall: T = 500 K
        % -----------------------------
        if j == 1
            Dsb = k*As/(0.5*dy);
            aB = Dsb;

            bc(i,j) = bc(i,j) + aB*T_bottom;
            ac_extra = ac_extra + aB;

            as(i,j) = 0;
        end

        % -----------------------------
        % North/top boundary
        % top wall: T = 300 K
        % top-right inlet: T = 300 K
        % -----------------------------
        if j == Ny
            Dnb = k*An/(0.5*dy);

            if i >= Nx - velocity_entry + 1
                % inlet: prescribed incoming temperature
                aB = Dnb + max(-FnT,0);
                T_b = T_inlet;
            else
                % top wall: prescribed T = 300
                aB = Dnb;
                T_b = T_top;
            end

            bc(i,j) = bc(i,j) + aB*T_b;
            ac_extra = ac_extra + aB;

            an(i,j) = 0;
        end

        ac(i,j) = ae(i,j) + aw(i,j) + an(i,j) + as(i,j) + ac_extra;

    end
end

end