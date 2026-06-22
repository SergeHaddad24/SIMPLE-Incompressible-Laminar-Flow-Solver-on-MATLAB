function [T_X,T_Y,T_up_X,T_up_Y,E_X,E_Y,E_up_X,E_up_Y,e_X,e_Y,e_up_X,e_up_Y,ce_distance_X,mag_ce,mag_ce_up,ce_distance_Y,ce_distance_up_X,ce_distance_up_Y] = ce_distance(m,n,centroid_X,centroid_Y,f_centroid_east_x,f_centroid_east_y,f_centroid_west_x,f_centroid_west_y,f_centroid_top_x,f_centroid_top_y,f_centroid_bottom_x,f_centroid_bottom_y,Swx,Swy,Sex,Sey,Snx,Sny,Ssx,Ssy)
ce_distance_X = zeros(m,n-1);
ce_distance_Y = zeros(m,n-1);
ce_distance_up_X= zeros(m-1,n);
ce_distance_up_Y= zeros(m-1,n);
e_X = zeros(m,n-1);
e_Y = zeros(m,n-1);
e_up_X= zeros(m-1,n);
e_up_Y = zeros(m-1,n);
mag_ce= zeros(m,n-1);
mag_ce_up= zeros(m-1,n);
E_X=zeros(m,n-1);
E_Y =zeros(m,n-1);
E_up_X=zeros(m-1,n);
E_up_Y=zeros(m-1,n);
T_X=zeros(m,n-1);
T_Y =zeros(m,n-1);
T_up_X=zeros(m-1,n);
T_up_Y=zeros(m-1,n);
for j =1:n-1
    for i=2:m-1
        ce_distance_X(i,j) = centroid_X(i,j)-centroid_X(i-1,j);
        ce_distance_Y(i,j) = centroid_Y(i,j)-centroid_Y(i-1,j);
        mag_ce(i,j)= sqrt(ce_distance_Y(i,j)^2+ce_distance_X(i,j)^2);
        e_X(i,j)=ce_distance_X(i,j)/mag_ce(i,j);
        e_Y(i,j)=ce_distance_Y(i,j)/mag_ce(i,j);
        E_X(i,j)=((Swx(i,j)^2 + Swy(i,j)^2)/(e_X(i,j)*Swx(i,j)+e_Y(i,j)*Swy(i,j)))*e_X(i,j);
        E_Y(i,j)=((Swx(i,j)^2 + Swy(i,j)^2)/(e_X(i,j)*Swx(i,j)+e_Y(i,j)*Swy(i,j)))*e_Y(i,j);
        T_X(i,j)=Swx(i,j)-E_X(i,j);
        T_Y(i,j)=Swy(i,j)-E_Y(i,j);
    end
    for i=1
        ce_distance_X(i,j) = centroid_X(i,j)-f_centroid_west_x(i,j);
        ce_distance_Y(i,j) = centroid_Y(i,j)-f_centroid_west_y(i,j);
        mag_ce(i,j)= sqrt(ce_distance_Y(i,j)^2+ce_distance_X(i,j)^2);
        e_X(i,j)=ce_distance_X(i,j)/mag_ce(i,j);
        e_Y(i,j)=ce_distance_Y(i,j)/mag_ce(i,j);
        E_X(i,j)=((Swx(i,j)^2 + Swy(i,j)^2)/(e_X(i,j)*Swx(i,j)+e_Y(i,j)*Swy(i,j)))*e_X(i,j);
        E_Y(i,j)=((Swx(i,j)^2 + Swy(i,j)^2)/(e_X(i,j)*Swx(i,j)+e_Y(i,j)*Swy(i,j)))*e_Y(i,j);
        T_X(i,j)=Swx(i,j)-E_X(i,j);
        T_Y(i,j)=Swy(i,j)-E_Y(i,j);
    end
    for i=m
        ce_distance_X(i,j)=f_centroid_east_x(i-1,j)-centroid_X(i-1,j);
        ce_distance_Y(i,j)=f_centroid_east_y(i-1,j)-centroid_Y(i-1,j);
        mag_ce(i,j)= sqrt(ce_distance_Y(i,j)^2+ce_distance_X(i,j)^2);
        e_X(i,j)=ce_distance_X(i,j)/mag_ce(i,j);
        e_Y(i,j)=ce_distance_Y(i,j)/mag_ce(i,j);
        E_X(i,j)=((Sex(i-1,j)^2 + Sey(i-1,j)^2)/(e_X(i,j)*Sex(i-1,j)+e_Y(i,j)*Sey(i-1,j)))*e_X(i,j);
        E_Y(i,j)=((Sex(i-1,j)^2 + Sey(i-1,j)^2)/(e_X(i,j)*Sex(i-1,j)+e_Y(i,j)*Sey(i-1,j)))*e_Y(i,j);
        T_X(i,j)=Sex(i-1,j)-E_X(i,j);
        T_Y(i,j)=Sey(i-1,j)-E_Y(i,j);
    end 
end
for i =1:m-1
    for j=2:n-1
        ce_distance_up_X(i,j) = centroid_X(i,j)-centroid_X(i,j-1);
        ce_distance_up_Y(i,j) = centroid_Y(i,j)-centroid_Y(i,j-1);
        mag_ce_up(i,j)= sqrt(ce_distance_up_Y(i,j)^2+ce_distance_up_X(i,j)^2);
        e_up_X(i,j)=ce_distance_up_X(i,j)/mag_ce_up(i,j);
        e_up_Y(i,j)=ce_distance_up_Y(i,j)/mag_ce_up(i,j);
        E_up_X(i,j)=((Ssx(i,j)^2 + Ssy(i,j)^2)/(e_up_X(i,j)*Ssx(i,j)+e_up_Y(i,j)*Ssy(i,j)))*e_up_X(i,j);
        E_up_Y(i,j)=((Ssx(i,j)^2 + Ssy(i,j)^2)/(e_up_X(i,j)*Ssx(i,j)+e_up_Y(i,j)*Ssy(i,j)))*e_up_Y(i,j);
        T_up_X(i,j)=Ssx(i,j)-E_up_X(i,j);
        T_up_Y(i,j)=Ssy(i,j)-E_up_Y(i,j);
    end
    for j=1
        ce_distance_up_X(i,j) = centroid_X(i,j)-f_centroid_bottom_x(i,j);
        ce_distance_up_Y(i,j) = centroid_Y(i,j)-f_centroid_bottom_y(i,j);
        mag_ce_up(i,j)= sqrt(ce_distance_up_Y(i,j)^2+ce_distance_up_X(i,j)^2);
        e_up_X(i,j)=ce_distance_up_X(i,j)/mag_ce_up(i,j);
        e_up_Y(i,j)=ce_distance_up_Y(i,j)/mag_ce_up(i,j);
        E_up_X(i,j)=((Ssx(i,j)^2 + Ssy(i,j)^2)/(e_up_X(i,j)*Ssx(i,j)+e_up_Y(i,j)*Ssy(i,j)))*e_up_X(i,j);
        E_up_Y(i,j)=((Ssx(i,j)^2 + Ssy(i,j)^2)/(e_up_X(i,j)*Ssx(i,j)+e_up_Y(i,j)*Ssy(i,j)))*e_up_Y(i,j);
        T_up_X(i,j)=Ssx(i,j)-E_up_X(i,j);
        T_up_Y(i,j)=Ssy(i,j)-E_up_Y(i,j);
    end
    for j=n
        ce_distance_up_X(i,j)=f_centroid_top_x(i,j-1)-centroid_X(i,j-1);
        ce_distance_up_Y(i,j)=f_centroid_top_y(i,j-1)-centroid_Y(i,j-1);
        mag_ce_up(i,j)= sqrt(ce_distance_up_Y(i,j)^2+ce_distance_up_X(i,j)^2);
        e_up_X(i,j)=ce_distance_up_X(i,j)/mag_ce_up(i,j);
        e_up_Y(i,j)=ce_distance_up_Y(i,j)/mag_ce_up(i,j);
        E_up_X(i,j)=((Snx(i,j-1)^2 + Sny(i,j-1)^2)/(e_up_X(i,j)*Snx(i,j-1)+e_up_Y(i,j)*Sny(i,j-1)))*e_up_X(i,j);
        E_up_Y(i,j)=((Snx(i,j-1)^2 + Sny(i,j-1)^2)/(e_up_X(i,j)*Snx(i,j-1)+e_up_Y(i,j)*Sny(i,j-1)))*e_up_Y(i,j);
        T_up_X(i,j)=Snx(i,j-1)-E_up_X(i,j);
        T_up_Y(i,j)=Sny(i,j-1)-E_up_Y(i,j);
    end 
end
end


       
