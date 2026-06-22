function [u_0, v_0, p_0, T_0 ] = initial_conditions(m,n)
u_0 = zeros(m-1,n-1);
v_0 = zeros(m-1,n-1);
p_0 = (10^5)*ones(m-1,n-1);
T_0 = 300*ones(m-1,n-1);
delta_x = 2/(m-1);
velocity_entry = round(0.2/delta_x);

for j=n-1
    for i=(m-1)-velocity_entry+1:m-1
        v_0(i,j)=-5;
        u_0(i,j)=0;
    end
end
for j=1
    for i=1:m-1
        T_0(i,j)=500;
    end
end
for j=n-1
    for i=1:m-1
        T_0(i,j)=300;
    end
end

end
