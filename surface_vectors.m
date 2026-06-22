function [Sex,Sey, Swx,Swy,Snx,Sny,Ssx,Ssy,nnx,nny,nex,ney] = surface_vectors(X,Y,m,n)
Sex=zeros(m-1,n-1);
Sey=zeros(m-1,n-1);
Swx=zeros(m-1,n-1);
Swy=zeros(m-1,n-1);
Snx=zeros(m-1,n-1);
Sny=zeros(m-1,n-1);
Ssx=zeros(m-1,n-1);
Ssy=zeros(m-1,n-1);
nwx=zeros(m-1,n-1);
nwy=zeros(m-1,n-1);
nnx=zeros(m-1,n-1);
nny=zeros(m-1,n-1);
nsx=zeros(m-1,n-1);
nsy=zeros(m-1,n-1);
nex=zeros(m-1,n-1);
ney=zeros(m-1,n-1);
for i=1:m-1
    for j=1:n-1
        Sex(i,j) = Y(i+1,j+1)-Y(i+1,j);
        Sey(i,j) = -(X(i+1,j+1)-X(i+1,j));
        nex(i,j) =Sex(i,j)/hypot(Sex(i,j),Sey(i,j));
        ney(i,j) =Sey(i,j)/hypot(Sex(i,j),Sey(i,j));
    end
end
for i=1:m-1
    for j=1:n-1
        Snx(i,j) = Y(i,j+1)-Y(i+1,j+1);
        Sny(i,j) = -(X(i,j+1)-X(i+1,j+1));
        nnx(i,j) =Snx(i,j)/hypot(Sex(i,j),Sey(i,j));
        nny(i,j) =Sny(i,j)/hypot(Sex(i,j),Sey(i,j));
    end
end
for i=1:m-1
    for j=1:n-1
        Swx(i,j) = Y(i,j)-Y(i,j+1);
        Swy(i,j) = -(X(i,j)-X(i,j+1));
        nwx(i,j) =Swx(i,j)/hypot(Sex(i,j),Sey(i,j));
        nwy(i,j) =Swy(i,j)/hypot(Sex(i,j),Sey(i,j));
    end
end
for i=1:m-1
    for j=1:n-1
        Ssx(i,j) = Y(i+1,j)-Y(i,j);
        Ssy(i,j) = -(X(i+1,j)-X(i,j));
        nsx(i,j) =Ssx(i,j)/hypot(Sex(i,j),Sey(i,j));
        nsy(i,j) =Ssy(i,j)/hypot(Sex(i,j),Sey(i,j));
    end
end
end



