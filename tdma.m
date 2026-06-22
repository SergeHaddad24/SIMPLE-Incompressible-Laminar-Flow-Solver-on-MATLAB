function [phi] = tdma(Nx, Ny, aC, aE, aW, aN, aS, bC, phi_old, urf)
    phi = phi_old;

    % Row sweep
    for j = 1:Ny
        a = zeros(Nx, 1);   % main
        b = zeros(Nx, 1);   % upper
        c = zeros(Nx, 1);   % lower
        d = zeros(Nx, 1);   % RHS

        for i = 1:Nx
            a(i) = aC(i,j) / urf;
            b(i) = aE(i,j);
            c(i) = aW(i,j);

            d(i) = bC(i,j);
            if j < Ny, d(i) = d(i) - aN(i,j) * phi(i,j+1); end
            if j > 1,  d(i) = d(i) - aS(i,j) * phi(i,j-1); end

            d(i) = d(i) + (1-urf)/urf * aC(i,j) * phi_old(i,j);
        end

        phi_line = tdma1D(a,b,c,d,Nx);
        phi(:,j) = phi_line;
    end

    % Column sweep
    for i = 1:Nx
        a = zeros(Ny, 1);   % main
        b = zeros(Ny, 1);   % upper
        c = zeros(Ny, 1);   % lower
        d = zeros(Ny, 1);   % RHS

        for j = 1:Ny
            a(j) = aC(i,j) / urf;
            b(j) = aN(i,j);
            c(j) = aS(i,j);

            d(j) = bC(i,j);
            if i < Nx, d(j) = d(j) - aE(i,j) * phi(i+1,j); end
            if i > 1,  d(j) = d(j) - aW(i,j) * phi(i-1,j); end

            d(j) = d(j) + (1-urf)/urf * aC(i,j) * phi_old(i,j);
        end

        phi_line = tdma1D(a,b,c,d,Ny);
        phi(i,:) = phi_line';
    end
end

function x = tdma1D(a,b,c,d,N)
    P = zeros(N, 1);
    Q = zeros(N, 1);
    x = zeros(N, 1);

    P(1) = -b(1) / a(1);
    Q(1) =  d(1) / a(1);

    for i = 2:N
        denom = a(i) + c(i) * P(i-1);
        P(i) = -b(i) / denom;
        Q(i) = (d(i) - c(i) * Q(i-1)) / denom;
    end

    x(N) = Q(N);
    for i = N-1:-1:1
        x(i) = P(i) * x(i+1) + Q(i);
    end
end