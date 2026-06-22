function p_prime = solve_pressure_correction(p_prime0,ae_pprime,aw_pprime,as_pprime,an_pprime,ac_pprime,bc_pprime,m,n,maxIter,tol)
p_prime = p_prime0;
urf_p = 0.5;   
dy = 1/(n-1);
% pressure_exit = round(0.1/dy);
for iterGS = 1:maxIter
    pold=p_prime;
    for i=1:m-1
        for j=1:n-1
           
            pE=0;pW=0;pS=0;pN=0;
            if i<m-1
                pE = p_prime(i+1,j);
            end
            if i>1
                pW=p_prime(i-1,j);
            end
            if j<n-1
                pN=p_prime(i,j+1);
            end
            if j>1 
                pS=p_prime(i,j-1);
            end
            p_star = (bc_pprime(i,j)+ae_pprime(i,j)*pE+aw_pprime(i,j)*pW+as_pprime(i,j)*pS+an_pprime(i,j)*pN)/ac_pprime(i,j);
            p_prime(i,j) = urf_p*p_star + (1-urf_p)*pold(i,j);
        end
    end
   
    res = max(abs(p_prime(:)-pold(:)));
    if res<tol
        break
    end
end
end
