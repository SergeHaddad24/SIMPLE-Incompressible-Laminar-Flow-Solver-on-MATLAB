function [phie,phiw,phis,phin] = geoInterp(m,n,gw,gs,gn,ge,phiC)
phie = zeros(m-1, n-1);
phiw = zeros(m-1, n-1);
phis = zeros(m-1, n-1);
phin = zeros(m-1, n-1);
for i=2:m-2
    for j=2:n-2
        phie(i,j)=ge(i,j)*phiC(i,j)+(1-ge(i,j))*phiC(i+1,j);
        phin(i,j)=gn(i,j)*phiC(i,j)+(1-gn(i,j))*phiC(i,j+1);
        phiw(i,j) = gw(i,j)*phiC(i,j) + (1-gw(i,j))*phiC(i-1,j);
        phis(i,j) = gs(i,j)*phiC(i,j) + (1-gs(i,j))*phiC(i,j-1);
    end
    for j=1
        phie(i,j)=ge(i,j)*phiC(i,j)+(1-ge(i,j))*phiC(i+1,j);
        phin(i,j)=gn(i,j)*phiC(i,j)+(1-gn(i,j))*phiC(i,j+1);
        phiw(i,j) = gw(i,j)*phiC(i,j) + (1-gw(i,j))*phiC(i-1,j);
        phis(i,j) = 320;
    end
    for j=n-1
        phie(i,j)=ge(i,j)*phiC(i,j)+(1-ge(i,j))*phiC(i+1,j);
        phis(i,j) = gs(i,j)*phiC(i,j) + (1-gs(i,j))*phiC(i,j-1);
        phiw(i,j) = gw(i,j)*phiC(i,j) + (1-gw(i,j))*phiC(i-1,j);
        phin(i,j) = phiC(i,j);
    end

end
for i=1
    for j=2:n-2
        phie(i,j)=ge(i,j)*phiC(i,j)+(1-ge(i,j))*phiC(i+1,j);
        phin(i,j)=gn(i,j)*phiC(i,j)+(1-gn(i,j))*phiC(i,j+1);
        phis(i,j) = gs(i,j)*phiC(i,j) + (1-gs(i,j))*phiC(i,j-1);
        phiw(i,j) = 400;
    end
    for j=1
        phie(i,j)=ge(i,j)*phiC(i,j)+(1-ge(i,j))*phiC(i+1,j);
        phin(i,j)=gn(i,j)*phiC(i,j)+(1-gn(i,j))*phiC(i,j+1);
        phis(i,j)=320;
        phiw(i,j)=400;
    end
    for j=n-1
        phie(i,j)=ge(i,j)*phiC(i,j)+(1-ge(i,j))*phiC(i+1,j);
        phis(i,j) = gs(i,j)*phiC(i,j) + (1-gs(i,j))*phiC(i,j-1);
        phiw(i,j)=400;
        phin(i,j)=phiC(i,j);
    end
end
for i =m-1
    for j=2:n-2
        phin(i,j)=gn(i,j)*phiC(i,j)+(1-gn(i,j))*phiC(i,j+1);
        phiw(i,j) = gw(i,j)*phiC(i,j) + (1-gw(i,j))*phiC(i-1,j);
        phis(i,j) = gs(i,j)*phiC(i,j) + (1-gs(i,j))*phiC(i,j-1);
        phie(i,j)=phiC(i,j);
    end
    for j=1
        phin(i,j)=gn(i,j)*phiC(i,j)+(1-gn(i,j))*phiC(i,j+1);
        phiw(i,j) = gw(i,j)*phiC(i,j) + (1-gw(i,j))*phiC(i-1,j);
        phis(i,j)=320;
        phie(i,j)=phiC(i,j);
    end
    for j=n-1
        phis(i,j) = gs(i,j)*phiC(i,j) + (1-gs(i,j))*phiC(i,j-1);
        phiw(i,j) = gw(i,j)*phiC(i,j) + (1-gw(i,j))*phiC(i-1,j);
        phie(i,j) = phiC(i,j);
        phin(i,j) = phiC(i,j);
    end
end
end


