function [S_c]= S(T,centroid_X,centroid_Y,m,n)
S_c=zeros(m-1,n-1);
for i =1:m-1
    for j=1:n-1
        S_c(i,j)= (2*centroid_X(i,j)-0.2*centroid_Y(i,j))/400;
    end
end
