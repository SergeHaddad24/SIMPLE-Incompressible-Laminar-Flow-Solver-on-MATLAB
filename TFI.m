function [X,Y]=TFI(m,n)
% number of discretizations along xi and eta axis



% discretize along xi and eta axis
xi = linspace(0.,1,m) ;
eta = linspace(0.,1.,n) ;

% Initialize matrices in x and y axis
X = zeros(m,n) ;
Y = zeros(m,n) ;


for i = 1:m
    Xi = xi(i) ;
    for j = 1:n
        Eta = eta(j) ;
        
        % Transfinite Interpolation 
        XY = (1-Eta)*Xb(Xi)+Eta*Xt(Xi)+(1-Xi)*Xl(Eta)+Xi*Xr(Eta)......
            -(Xi*Eta*Xt(1)+Xi*(1-Eta)*Xb(1)+Eta*(1-Xi)*Xt(0)+(1-Xi)*(1-Eta)*Xb(0)) ;
    
        X(i,j) = XY(1) ;
        Y(i,j) = XY(2) ;
        
    end
end

