function[gamma_C]=Gamma(centroid_X,centroid_Y,T,m,n)
gamma_C =zeros(m-1,n-1);
for i=1:m-1
    for j=1:n-1
        gamma_C(i,j)=(centroid_X(i,j)^2+exp(0.1*centroid_Y(i,j)))/400;
    end
end
