function [f_centroid_east_x,f_centroid_east_y,f_centroid_west_x,f_centroid_west_y,f_centroid_top_x,f_centroid_top_y,f_centroid_bottom_x,f_centroid_bottom_y] = face_centroid(X,Y,m,n)
f_centroid_east_x=zeros(m-1,n-1);
f_centroid_west_x=zeros(m-1,n-1);
f_centroid_top_x=zeros(m-1,n-1);
f_centroid_bottom_x=zeros(m-1,n-1);
f_centroid_east_y=zeros(m-1,n-1);
f_centroid_west_y=zeros(m-1,n-1);
f_centroid_top_y=zeros(m-1,n-1);
f_centroid_bottom_y=zeros(m-1,n-1);
for i=1:m-1
    for j=1:n-1
        f_centroid_east_x(i,j) =(X(i+1,j)+X(i+1,j+1))/2;
        f_centroid_west_x(i,j) = (X(i,j)+X(i,j+1))/2;
        f_centroid_top_x(i,j) = (X(i,j+1)+X(i+1,j+1))/2;
        f_centroid_bottom_x(i,j) = (X(i,j)+ X(i+1,j))/2;
        f_centroid_east_y(i,j) =(Y(i+1,j)+Y(i+1,j+1))/2;
        f_centroid_west_y(i,j) = (Y(i,j)+Y(i,j+1))/2;
        f_centroid_top_y(i,j) = (Y(i,j+1)+Y(i+1,j+1))/2;
        f_centroid_bottom_y(i,j) = (Y(i,j)+ Y(i+1,j))/2;
    end
end