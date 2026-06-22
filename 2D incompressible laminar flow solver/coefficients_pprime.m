function [ae_pprime,aw_pprime,as_pprime,an_pprime,ac_pprime,bc_pprime] = coefficients_pprime(Fe,Fw,Fn,Fs,ac_u,ac_v,rho,area,m,n,Sex,Sey,Swx,Swy,Snx,Sny,Ssx,Ssy,mag_ce,mag_ce_up)
ae_pprime=zeros(m-1,n-1);
aw_pprime=zeros(m-1,n-1);
as_pprime=zeros(m-1,n-1);
an_pprime=zeros(m-1,n-1);
ac_pprime=zeros(m-1,n-1);
bc_pprime=zeros(m-1,n-1);
eps0=1e-14;
dy = 1/(n-1);
pressure_exit = round(0.1/dy);

for i=1:m-1
    for j=1:n-1
        aW_boundary=0;
        Ae= hypot(Sex(i,j),Sey(i,j));
        Aw= hypot(Swx(i,j),Swy(i,j));
        As= hypot(Ssx(i,j),Ssy(i,j));
        An= hypot(Snx(i,j),Sny(i,j));
        if i<m-1
            dP=area(i,j)/max(abs(ac_u(i,j)),eps0);
            dE=area(i+1,j)/max(abs(ac_u(i+1,j)),eps0);
            d_face= 0.5*(dP+dE);
            ae_pprime(i,j)=rho*d_face*Ae/mag_ce(i+1,j);
        else
            ae_pprime(i,j)=0;
        end
        if i>1
            dP=area(i,j)/max(abs(ac_u(i,j)),eps0);
            dW=area(i-1,j)/max(abs(ac_u(i-1,j)),eps0);
            d_face= 0.5*(dP+dW);
            aw_pprime(i,j)=rho*d_face*Aw/mag_ce(i,j);

        else
            aw_pprime(i,j)=0;
            if j <= pressure_exit
                % west pressure outlet, p'_boundary = 0
                dP = area(i,j)/max(abs(ac_u(i,j)),eps0);
                aW_boundary = rho*dP*Aw/mag_ce(i,j);
            else
                aW_boundary = 0;
            end
        end
        if j<n-1
            dP=area(i,j)/max(abs(ac_v(i,j)),eps0);
            dN=area(i,j+1)/max(abs(ac_v(i,j+1)),eps0);
            d_face= 0.5*(dP+dN);
            an_pprime(i,j)=rho*d_face*An/mag_ce_up(i,j+1);
        else
            an_pprime(i,j)=0;
            
        end
        
        if j>1
            dP=area(i,j)/max(abs(ac_v(i,j)),eps0);
            dS=area(i,j-1)/max(abs(ac_v(i,j-1)),eps0);
            d_face= 0.5*(dP+dS);
            as_pprime(i,j)=rho*d_face*As/mag_ce_up(i,j);
        else
            as_pprime(i,j)=0;
        end
        if i == 1 && j <= pressure_exit
            ac_pprime(i,j) = ae_pprime(i,j)+aw_pprime(i,j)+an_pprime(i,j)+as_pprime(i,j) + aW_boundary;
        else
            ac_pprime(i,j) = ae_pprime(i,j)+aw_pprime(i,j)+an_pprime(i,j)+as_pprime(i,j);
        end
        bc_pprime(i,j)=-(Fe(i,j)+Fw(i,j)+Fs(i,j)+Fn(i,j));
    end
end

end