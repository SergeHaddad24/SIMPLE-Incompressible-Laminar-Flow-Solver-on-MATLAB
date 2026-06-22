function settings = getinfo()

    fprintf('--- Solver Settings ---\n');

    settings.Nx = input('Enter number of cells in x-direction: ');
    settings.Ny = input('Enter number of cells in y-direction: ');

   

    settings.urf = input('Enter under-relaxation factor (e.g. 0.8): ');
    settings.tol = input('Enter convergence tolerance (e.g. 1e-6): ');
    settings.maxIter = input('Enter maximum number of iterations: ');

    
end