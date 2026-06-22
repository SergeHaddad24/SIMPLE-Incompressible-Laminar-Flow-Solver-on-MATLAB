function phi_f = smart_face(phiU,phiC,phiD)

eps0 = 1e-14;
den = phiC - phiD;

if abs(den) < eps0
    phi_f = phiC;
    return
end

phihatU = (phiU - phiD)/den;

if phihatU <= 0
    phihatF = 0;
elseif phihatU <= 1/6
    phihatF = 3*phihatU;
elseif phihatU <= 5/6
    phihatF = (3/8)*(2*phihatU + 1);
elseif phihatU <= 1
    phihatF = 1;
else
    phihatF = phihatU;
end

phi_f = phihatF*(phiC - phiD) + phiD;
end