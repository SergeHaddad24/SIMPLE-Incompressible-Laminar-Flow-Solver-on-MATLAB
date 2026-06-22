function [area] = get_area(X,Y,m,n)
area=zeros(m-1,n-1);
area1=zeros(m-1,n-1);
area2=zeros(m-1,n-1);
for i=1:m-1
    for j=1:n-1
        area1(i,j)=1/2*((X(i+1,j)-X(i,j))*(Y(i+1,j+1)-Y(i,j))-(Y(i+1,j)-Y(i,j))*(X(i+1,j+1)-X(i,j)));
        area2(i,j)=1/2*((X(i+1,j+1)-X(i,j))*(Y(i,j+1)-Y(i,j))-(Y(i+1,j+1)-Y(i,j))*(X(i,j+1)-X(i,j)));
        area(i,j) = area1(i,j) + area2(i,j);
    end
end
end
