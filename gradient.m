function [dphidx,dphidy]=gradient(m,n,phie,phis,phiw,phin,area,Ssx,Ssy,Sex,Sey,Snx,Sny,Swx,Swy)
dphidx = zeros(m-1,n-1);
dphidy = zeros(m-1,n-1);
sumX = zeros(m-1,n-1);
sumY= zeros(m-1,n-1);
for i=1:m-1
    for j=1:n-1
        sumX(i,j)=phie(i,j)*Sex(i,j)+phiw(i,j)*Swx(i,j)+phin(i,j)*Snx(i,j)+phis(i,j)*Ssx(i,j);
        sumY(i,j)=phin(i,j)*Sny(i,j)+phiw(i,j)*Swy(i,j)+phie(i,j)*Sey(i,j)+phis(i,j)*Ssy(i,j);
        dphidx(i,j) = sumX(i,j) / area(i,j);
        dphidy(i,j) = sumY(i,j) / area(i,j);
    end
end