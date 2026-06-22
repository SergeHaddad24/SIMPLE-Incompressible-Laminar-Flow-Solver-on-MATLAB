function phi = solve_GS_scalar(phi,ae,aw,as,an,ac,bc,m,n,urf,maxIter,tol)

for iter = 1:maxIter
    phi_old = phi;

    for i = 1:m-1
        for j = 1:n-1

            phiE = 0; phiW = 0; phiN = 0; phiS = 0;

            if i < m-1, phiE = phi(i+1,j); end
            if i > 1,   phiW = phi(i-1,j); end
            if j < n-1, phiN = phi(i,j+1); end
            if j > 1,   phiS = phi(i,j-1); end

            phi_star = (bc(i,j) + ae(i,j)*phiE + aw(i,j)*phiW + ...
                        an(i,j)*phiN + as(i,j)*phiS)/ac(i,j);

            phi(i,j) = urf*phi_star + (1-urf)*phi_old(i,j);
        end
    end

    res = max(abs(phi(:)-phi_old(:)));

    if mod(iter,100) == 0
        fprintf('T GS iter %d, residual = %.6e\n',iter,res);
    end

    if res < tol
        fprintf('T converged in %d iterations, residual = %.6e\n',iter,res);
        break
    end
end

end