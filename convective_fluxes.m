function [Fe, Fw, Fn, Fs] = convective_fluxes(rho,m,n,u,v,Sex,Sey,Swx,Swy,Snx,Sny,Ssx,Ssy)
Fe =zeros(m-1,n-1);
Fw =zeros(m-1,n-1);
Fs=zeros(m-1,n-1);
Fn=zeros(m-1,n-1);
uf=0;
vf=0;
velocity_entry = round(0.2/(2/(m-1)));
for i=1:m-1
    for j=1:n-1
                
        
        
        if i<m-1
            uf = 0.5*(u(i,j) + u(i+1,j));
            vf = 0.5*(v(i,j) + v(i+1,j));
        else 
            uf=0;
            vf=0;
        end
        Fe(i,j)=rho*(uf*Sex(i,j)+vf*Sey(i,j));

       
        if i>1
            uf = 0.5*(u(i-1,j) + u(i,j));
            vf = 0.5*(v(i-1,j) + v(i,j));
        else 
            uf=0;
            vf=0;
        end
         Fw(i,j) = rho*(uf*Swx(i,j) + vf*Swy(i,j));
        if j<n-1
            uf = 0.5*(u(i,j) + u(i,j+1));
            vf = 0.5*(v(i,j) + v(i,j+1));
        else 
            
            if i>= (m-1)-velocity_entry+1
                uf=0;
                vf=-5;
            else
                uf=0;
                vf=0;
            end


            
        end
        Fn(i,j) = rho*(uf*Snx(i,j) + vf*Sny(i,j));
        if j>1
            uf = 0.5*(u(i,j-1) + u(i,j));
            vf = 0.5*(v(i,j-1) + v(i,j));
        else
            uf =0;
            vf=0;
        end
        Fs(i,j) = rho*(uf*Ssx(i,j) + vf*Ssy(i,j));
    end
end
end
