function phi = SMART(phi, ae, aw, as, an, ac, bc_base, Fe, Fw, Fn, Fs, m, n, urf, maxIter, maxOuter, tol)


phi_prev = phi;                       
Sp = (1 - urf) .* ac .* phi_prev;     
beta = 0.5;                          

% ---------- upwind warm-up (pure Gauss-Seidel) ----------

for iterGS = 1:maxIter
    phioldGS = phi;
    for i = 1:m-1
        for j = 1:n-1
            phiE = 0.0; phiW = 0.0; phiN = 0.0; phiS = 0.0;
            if i < m-1, phiE = phi(i+1,j); end
            if i > 1,   phiW = phi(i-1,j); end
            if j < n-1, phiN = phi(i,j+1); end
            if j > 1,   phiS = phi(i,j-1); end
            phi(i,j) = (bc_base(i,j) + Sp(i,j) + ae(i,j)*phiE + aw(i,j)*phiW + ...
                        an(i,j)*phiN + as(i,j)*phiS) / ac(i,j);
        end
    end
    if max(abs(phi(:) - phioldGS(:))) < tol, break; end
end

% ---------- SMART deferred-correction outer loop ----------
for iterOuter = 1:max(maxOuter,1)
    phiOldOuter = phi;

    
    bdc = deferred_correction_smart(phiOldOuter, m, n, Fe, Fw, Fn, Fs);
    bc  = bc_base + bdc + Sp;         

    phiSolve = phi;
    for iterGS = 1:maxIter
        phioldGS = phiSolve;
        for i = 1:m-1
            for j = 1:n-1
                phiE = 0.0; phiW = 0.0; phiN = 0.0; phiS = 0.0;
                if i < m-1, phiE = phiSolve(i+1,j); end
                if i > 1,   phiW = phiSolve(i-1,j); end
                if j < n-1, phiN = phiSolve(i,j+1); end
                if j > 1,   phiS = phiSolve(i,j-1); end
                phiSolve(i,j) = (bc(i,j) + ae(i,j)*phiE + aw(i,j)*phiW + ...
                                 an(i,j)*phiN + as(i,j)*phiS) / ac(i,j);
            end
        end
        if max(abs(phiSolve(:) - phioldGS(:))) < tol, break; end
    end

  
    phi = phiOldOuter + beta * (phiSolve - phiOldOuter);

    resOuter = max(abs(phi(:) - phiOldOuter(:)));
    fprintf('SMART outer iter = %d, res = %.6e\n', iterOuter, resOuter);
    if resOuter < tol, break; end
end
end
