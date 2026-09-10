function phi_f = smart_face(phiU,phiC,phiD)

eps0 = 1e-14;
den = phiD - phiU;

if abs(den) < eps0
    phi_f = phiC;
    return
end

phihatC = (phiC - phiU)/den;

if phihatC <= 0 || phihatC >= 1
    phihatF = phihatC;
elseif phihatC < 1/6
    phihatF = 3*phihatC;
elseif phihatC < 5/6
    phihatF = (3/8)*(2*phihatC + 1);
else
    phihatF = 1;
end

phi_f = phihatF*den + phiU;
end
