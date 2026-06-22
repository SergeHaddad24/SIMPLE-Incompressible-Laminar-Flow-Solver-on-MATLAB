function [ge,gw,gs,gn] =get_interpolationFactors(m,n,area)
ge=zeros(m-1,n-1);
gw=zeros(m-1,n-1);
gs=zeros(m-1,n-1);
gn=zeros(m-1,n-1);
for i =2:m-2
    for j=2:n-2
        ge(i,j)=area(i,j)/(area(i,j)+area(i+1,j));
        gw(i,j)=area(i,j)/(area(i,j)+area(i-1,j));
        gs(i,j)=area(i,j)/(area(i,j)+area(i,j-1));
        gn(i,j)=area(i,j)/(area(i,j)+area(i,j+1));
    end
    for j=1
        ge(i,j) = area(i,j) / (area(i,j) + area(i+1,j));
        gw(i,j) = area(i,j) / (area(i,j) + area(i-1,j));
        gn(i,j) = area(i,j) / (area(i,j) + area(i,j+1));
        gs(i,j)=1;
    end
    for j=n-1
        ge(i,j)=area(i,j)/(area(i,j)+area(i+1,j));
        gw(i,j)=area(i,j)/(area(i,j)+area(i-1,j));
        gs(i,j)=area(i,j)/(area(i,j)+area(i,j-1));
        gn(i,j)=1;
    end
end
for i=1
     for j=2:n-2
        ge(i,j)=area(i,j)/(area(i,j)+area(i+1,j));
        gs(i,j)=area(i,j)/(area(i,j)+area(i,j-1));
        gn(i,j)=area(i,j)/(area(i,j)+area(i,j+1));
        gw(i,j)=1;
    end
    for j=1
        ge(i,1) = area(i,1) / (area(i,1) + area(i+1,1));
        gn(i,1) = area(i,1) / (area(i,1) + area(i,1+1));
        gs(i,j)=1;
        gw(i,j)=1;
    end
    for j=n-1
        ge(i,j)=area(i,j)/(area(i,j)+area(i+1,j));
        gs(i,j)=area(i,j)/(area(i,j)+area(i,j-1)); 
        gw(i,j)=1;
        gn(i,j)=1;
    end
end
for i=m-1
    for j=2:n-2
        gw(i,j)=area(i,j)/(area(i,j)+area(i-1,j));
        gs(i,j)=area(i,j)/(area(i,j)+area(i,j-1));
        gn(i,j)=area(i,j)/(area(i,j)+area(i,j+1));
        ge(i,j)=1;
    end
    for j=1
        gw(i,1) = area(i,j) / (area(i,j) + area(i-1,j));
        gn(i,1) = area(i,j) / (area(i,j) + area(i,j+1));
        ge(i,j)=1;
        gs(i,j)=1;
    end
    for j=n-1
        gw(i,j)=area(i,j)/(area(i,j)+area(i-1,j));
        gs(i,j)=area(i,j)/(area(i,j)+area(i,j-1));
        gn(i,j)=1;
        ge(i,j)=1;
    end
end