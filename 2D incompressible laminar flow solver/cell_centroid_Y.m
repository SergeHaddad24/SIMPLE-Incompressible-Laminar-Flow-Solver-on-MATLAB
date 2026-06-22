function centroid_Y= cell_centroid_Y(Y,m,n)
 % c_x = (X(x,y)+ X(x+1,y)+X(x,y+1)+X(x+1,y+1))/4;
% c_y = (Y(x,y)+Y(x+1,y)+Y(x,y+1)+Y(x+1,y+1))/4;
 centroid_Y = zeros(m-1,n-1);
 for i=1:m-1
     for j=1:n-1
       centroid_Y(i,j)= (Y(i,j)+Y(i+1,j)+Y(i,j+1)+Y(i+1,j+1))/4;
     end
end
