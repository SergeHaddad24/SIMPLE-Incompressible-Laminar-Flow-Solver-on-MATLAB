function [ae,aw,as,an,ac,bc_base] = coefficients_v(du_dy,dv_dy,mu,p,Sex,Sey,Ssx,Ssy,Snx,Sny,Swx,Swy,E_X,E_Y,mag_ce,mag_ce_up,E_up_X,E_up_Y,area,m,n,Fe,Fw,Fn,Fs)
ae=zeros(m-1,n-1);
aw=zeros(m-1,n-1);
as=zeros(m-1,n-1);
an=zeros(m-1,n-1);
ac=zeros(m-1,n-1);
bc_base=zeros(m-1,n-1);
% ce=zeros(m-1,n-1);
% cw=zeros(m-1,n-1);
% cs=zeros(m-1,n-1);
% cn=zeros(m-1,n-1);
% ae0=zeros(m-1,n-1);
% aw0=zeros(m-1,n-1);
% as0=zeros(m-1,n-1);
% an0=zeros(m-1,n-1);



delta_y =1/(n-1);
% gradTw_x=zeros(m-1,n-1);
% gradTw_y=zeros(m-1,n-1);
% gradTs_x=zeros(m-1,n-1);
% gradTs_y=zeros(m-1,n-1);
% gradTe_x=zeros(m-1,n-1);
% gradTe_y=zeros(m-1,n-1);
% gradTn_x=zeros(m-1,n-1);
% gradTn_y=zeros(m-1,n-1);
velocity_entry = round(0.2/(2/(m-1)));
pressure_exit=round(0.1/delta_y);
for i=1:m-1
    for j=1:n-1
        % ae0(i,j)=-gamma_e(i,j)*(hypot(E_X(i+1,j),E_Y(i+1,j))/mag_ce(i+1,j));
        % aw0(i,j)=-gamma_w(i,j)*(hypot(E_X(i,j),E_Y(i,j))/mag_ce(i,j));
        % as0(i,j)=-gamma_s(i,j)*(hypot(E_up_X(i,j),E_up_Y(i,j))/mag_ce_up(i,j));
        % an0(i,j)=-gamma_n(i,j)*(hypot(E_up_X(i,j+1),E_up_Y(i,j+1))/mag_ce_up(i,j+1));
        % ae(i,j) = ae0(i,j);
        % an(i,j) = an0(i,j);
        % as(i,j) = as0(i,j);
        % aw(i,j) = aw0(i,j);
        Ae = hypot(Sex(i,j),Sey(i,j));
        Aw = hypot(Swx(i,j),Swy(i,j));
        An = hypot(Snx(i,j),Sny(i,j));
        As = hypot(Ssx(i,j),Ssy(i,j));

        De = mu*Ae/mag_ce(i+1,j);
        Dw = mu*Aw/mag_ce(i,j);
        Dn = mu*An/mag_ce_up(i,j+1);
        Ds = mu*As/mag_ce_up(i,j);
        ae(i,j)= max(-Fe(i,j),0)+De;
        aw(i,j)= max(-Fw(i,j),0)+Dw;
        as(i,j)= max(-Fs(i,j),0)+Ds;
        an(i,j)= max(-Fn(i,j),0)+Dn;
       
        ac_extra=0;
       
        
        
        
        if j == 1
            gradP_y = (p(i,j+1)-p(i,j))/delta_y;
        elseif j == n-1
            gradP_y = (p(i,j)-p(i,j-1))/delta_y;
        else
            gradP_y = (p(i,j+1)-p(i,j-1))/(2*delta_y);
        end
        
        % if i<m-1, ce(i,j)=-gamma_e(i,j)*(dphidx_e(i,j)*T_X(i+1,j)+dphidy_e(i,j)*T_Y(i+1,j));end
        % if i>1, cw(i,j)=-gamma_w(i,j)*(dphidx_w(i,j)*T_X(i,j)+dphidy_w(i,j)*T_Y(i,j)); end
        % if j>1,cs(i,j)=-gamma_s(i,j)*(dphidx_s(i,j)*T_up_X(i,j)+dphidy_s(i,j)*T_up_Y(i,j));end
        % if j<n-1,cn(i,j)=-gamma_n(i,j)*(dphidx_n(i,j)*T_up_X(i,j+1)+dphidy_n(i,j)*T_up_Y(i,j+1));end
        % bc(i,j)=ce(i,j)+cn(i,j)+cw(i,j)+cs(i,j)+area(i,j)*S_c(i,j);
        
        if i < m-1
            du_dy_e = 0.5*(du_dy(i,j) + du_dy(i+1,j));
            dv_dy_e = 0.5*(dv_dy(i,j) + dv_dy(i+1,j));
            

        else
            du_dy_e = du_dy(i,j);
            dv_dy_e = dv_dy(i,j);
            
        end

        
        if i > 1
            du_dy_w = 0.5*(du_dy(i-1,j) + du_dy(i,j));
            dv_dy_w = 0.5*(dv_dy(i-1,j) + dv_dy(i,j));
            
        else
            du_dy_w = du_dy(i,j);
            dv_dy_w = dv_dy(i,j);
            
        end

        
        if j < n-1
            du_dy_n = 0.5*(du_dy(i,j) + du_dy(i,j+1));
            dv_dy_n = 0.5*(dv_dy(i,j) + dv_dy(i,j+1));
        else
            du_dy_n = du_dy(i,j);
            dv_dy_n = dv_dy(i,j);
        end

        
        if j > 1
            du_dy_s = 0.5*(du_dy(i,j-1) + du_dy(i,j));
            dv_dy_s = 0.5*(dv_dy(i,j-1) + dv_dy(i,j));
        else
            du_dy_s = du_dy(i,j);
            dv_dy_s = dv_dy(i,j);
        end
        transpose_stress_source = mu*(du_dy_e*Sex(i,j)+dv_dy_e*Sey(i,j))+mu*(du_dy_w*Swx(i,j)+dv_dy_w*Swy(i,j))+mu*(du_dy_s*Ssx(i,j)+dv_dy_s*Ssy(i,j))+mu*(du_dy_n*Snx(i,j)+dv_dy_n*Sny(i,j));
        bc_base(i,j) = transpose_stress_source - gradP_y*area(i,j);
        if i==1
          
           % gradTw_x(i,j) = ((400 - T(1,j))/mag_ce(1,j)) * e_X(1,j);
           % gradTw_y(i,j) = ((400 - T(1,j))/mag_ce(1,j)) * e_Y(1,j);
           % cw(i,j) = -gamma_w(i,j) * (gradTw_x(i,j)*T_X(i,j) + gradTw_y(i,j)*T_Y(i,j));
           % bc(i,j) = bc(i,j) -aw0(i,j)*400;
           % aw(i,j)=0;
            
           if j <= pressure_exit
               % pressure outlet: zero-gradient velocity, no Dirichlet wall contribution
               aw(i,j) = 0.0;
           else
               % west wall above outlet
               A_w = hypot(Swx(i,j),Swy(i,j));
               dx = 2/(m-1);

               Dwb = mu*A_w/(0.5*dx);

               v_b = 0;

               bc_base(i,j) = bc_base(i,j) + Dwb*v_b;
               ac_extra = ac_extra + Dwb;

               aw(i,j) = 0.0;
           end

            
            
        end
        if j==1
           % as(i,j)=-gamma_s(i,j)*(hypot(E_up_X(i,j),E_up_Y(i,j))/mag_ce_up(i,j));
           % gradTs_x(i,j) = ((320 - T(i,j))/mag_ce_up(i,j)) * e_up_X(i,j);
           % gradTs_y(i,j) = ((320 - T(i,j))/mag_ce_up(i,j)) * e_up_Y(i,j);
           % cs(i,j) = -gamma_s(i,j) * (gradTs_x(i,j)*T_up_X(i,j) + gradTs_y(i,j)*T_up_Y(i,j));
           % bc(i,j) = bc(i,j) -as0(i,j)*320;
           % as(i,j)=0;
           
               A_s = hypot(Ssx(i,j),Ssy(i,j));
               dy = 1/(n-1);

               Dsb = mu*A_s/(0.5*dy);

               v_b = 0;

               bc_base(i,j) = bc_base(i,j) + Dsb*v_b;
               ac_extra = ac_extra + Dsb;

               as(i,j) = 0.0;
           
           
        end
        
        if j==n-1
            % gradTn_x(i,j)=dphidx_n(i,j)-(dphidx_n(i,j)*nnx(i,j)^2+dphidy_n(i,j)*nnx(i,j)*nny(i,j));
            % gradTn_y(i,j)=dphidy_n(i,j)-(dphidx_n(i,j)*nnx(i,j)*nny(i,j)+dphidy_n(i,j)*nny(i,j)^2);
            % % cn(i,j)=-gamma_n(i,j)*(gradTn_x(i,j)*T_up_X(i,j+1)+gradTn_y(i,j)*T_up_Y(i,j+1));
            % cn(i,j)=0;
            % an(i,j)=0;
            % an0(i,j)=0;
            if i>= (m-1)-velocity_entry+1
                v_b=-5;
                A_n = hypot(Snx(i,j),Sny(i,j));
                dy = 1/(n-1);

                Dnb = mu*A_n/(0.5*dy);

                Fnb = Fn(i,j);          
                aB  = Dnb + max(-Fnb,0);

                bc_base(i,j) = bc_base(i,j) + aB*v_b;
                ac_extra = ac_extra + aB;
            else
                v_b=0;
                A_n = hypot(Snx(i,j),Sny(i,j));
                dy = 1/(n-1);

                Dnb = mu*A_n/(0.5*dy);

                aB = Dnb;

                bc_base(i,j) = bc_base(i,j) + aB*v_b;
                ac_extra = ac_extra + aB;
            end
            
            an(i,j)=0;
            
        end
        % ac(i,j)=-(ae0(i,j)+an0(i,j)+as0(i,j)+aw0(i,j));
        if i==m-1
            % gradTe_x(i,j)=dphidx_e(i,j)-(dphidx_e(i,j)*nex(i,j)^2+dphidy_e(i,j)*nex(i,j)*ney(i,j)+(-(15/gamma_e(i,j)*(T(i,j)-300)))*nex(i,j));
            % gradTe_y(i,j)=dphidx_e(i,j)-(dphidx_e(i,j)*nex(i,j)*ney(i,j)+dphidy_e(i,j)*ney(i,j)^2+(-(15/gamma_e(i,j)*(T(i,j)-300)))*ney(i,j));
            % ce(i,j)=-gamma_e(i,j)*(gradTe_x(i,j)*T_X(i+1,j)+gradTe_y(i,j)*T_Y(i+1,j));
            % ae(i,j)=0;
            % Aface_e(i,j)=hypot(Sex(i,j),Sey(i,j));
            % bc(i,j) = bc(i,j) + 15*300*Aface_e(i,j);
            % ac(i,j) = ac(i,j) +15*Aface_e(i,j);
            A_e = hypot(Sex(i,j),Sey(i,j));
            dx = 2/(m-1);

            Deb = mu*A_e/(0.5*dx);

            v_b = 0;

            bc_base(i,j) = bc_base(i,j) + Deb*v_b;
            ac_extra = ac_extra + Deb;

            ae(i,j) = 0;


        end
         

         
         
         
         ac(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+ac_extra;
        
    end
end
