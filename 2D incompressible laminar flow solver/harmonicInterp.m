function [gamma_e,gamma_s,gamma_n,gamma_w]=harmonicInterp(m,n,ge,gs,gn,gw,gamma_C)
gamma_e=zeros(m-1,n-1);
gamma_s = zeros(m-1, n-1);
gamma_n = zeros(m-1, n-1);
gamma_w = zeros(m-1, n-1);
for i=2:m-2
    for j=2:n-2
        gamma_e(i,j)=1/(ge(i,j)/gamma_C(i,j)+(1-ge(i,j))/gamma_C(i+1,j));
        gamma_n(i,j)=1/(gn(i,j)/gamma_C(i,j)+(1-gn(i,j))/gamma_C(i,j+1));
        gamma_w(i,j) = 1/(gw(i,j)/gamma_C(i,j) + (1-gw(i,j))/gamma_C(i-1,j));
        gamma_s(i,j) = 1/(gs(i,j)/gamma_C(i,j) + (1-gs(i,j))/gamma_C(i,j-1));
    end
    for j=1
        gamma_e(i,j)=1/(ge(i,j)/gamma_C(i,j)+(1-ge(i,j))/gamma_C(i+1,j));
        gamma_n(i,j)=1/(gn(i,j)/gamma_C(i,j)+(1-gn(i,j))/gamma_C(i,j+1));
        gamma_w(i,j) = 1/(gw(i,j)/gamma_C(i,j) + (1-gw(i,j))/gamma_C(i-1,j));
        gamma_s(i,j)= gamma_C(i,j);
    end
    for j=n-1
        gamma_e(i,j)=1/(ge(i,j)/gamma_C(i,j)+(1-ge(i,j))/gamma_C(i+1,j));
        gamma_w(i,j) = 1/(gw(i,j)/gamma_C(i,j) + (1-gw(i,j))/gamma_C(i-1,j));
        gamma_s(i,j) = 1/(gs(i,j)/gamma_C(i,j) + (1-gs(i,j))/gamma_C(i,j-1));
        gamma_n(i,j) = gamma_C(i,j);
    end

end
for i=1
    for j=2:n-2
        gamma_e(i,j)=1/(ge(i,j)/gamma_C(i,j)+(1-ge(i,j))/gamma_C(i+1,j));
        gamma_n(i,j)=1/(gn(i,j)/gamma_C(i,j)+(1-gn(i,j))/gamma_C(i,j+1));
        gamma_s(i,j) = 1/(gs(i,j)/gamma_C(i,j) + (1-gs(i,j))/gamma_C(i,j-1));
        gamma_w(i,j) =gamma_C(i,j);
    end
    for j=1
        gamma_e(i,j)=1/(ge(i,j)/gamma_C(i,j)+(1-ge(i,j))/gamma_C(i+1,j));
        gamma_n(i,j)=1/(gn(i,j)/gamma_C(i,j)+(1-gn(i,j))/gamma_C(i,j+1));
        gamma_s(i,j) = gamma_C(i,j);
        gamma_w(i,j)=gamma_C(i,j);
    end
    for j=n-1
        gamma_e(i,j)=1/(ge(i,j)/gamma_C(i,j)+(1-ge(i,j))/gamma_C(i+1,j));
        gamma_s(i,j) = 1/(gs(i,j)/gamma_C(i,j) + (1-gs(i,j))/gamma_C(i,j-1));
        gamma_w(i,j) = gamma_C(i,j);
        gamma_n(i,j) = gamma_C(i,j);
    end
end
for i =m-1
    for j=2:n-2
        gamma_n(i,j)=1/(gn(i,j)/gamma_C(i,j)+(1-gn(i,j))/gamma_C(i,j+1));
        gamma_w(i,j) = 1/(gw(i,j)/gamma_C(i,j) + (1-gw(i,j))/gamma_C(i-1,j));
        gamma_s(i,j) = 1/(gs(i,j)/gamma_C(i,j) + (1-gs(i,j))/gamma_C(i,j-1));
        gamma_e(i,j) = gamma_C(i,j);
    end
    for j=1
        gamma_n(i,j)=1/(gn(i,j)/gamma_C(i,j)+(1-gn(i,j))/gamma_C(i,j+1));
        gamma_w(i,j) = 1/(gw(i,j)/gamma_C(i,j) + (1-gw(i,j))/gamma_C(i-1,j));
        gamma_e(i,j) = gamma_C(i,j);
        gamma_s(i,j) = gamma_C(i,j);
    end
    for j=n-1
        gamma_w(i,j) = 1/(gw(i,j)/gamma_C(i,j) + (1-gw(i,j))/gamma_C(i-1,j));
        gamma_s(i,j) = 1/(gs(i,j)/gamma_C(i,j) + (1-gs(i,j))/gamma_C(i,j-1));
        gamma_e(i,j) = gamma_C(i,j);
        gamma_n(i,j) = gamma_C(i,j);
    end
end
end