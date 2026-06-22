function bdc = deferred_correction_smart(phi,m,n,Fe,Fw,Fn,Fs)

bdc = zeros(m-1,n-1);

for i = 1:m-1
    for j = 1:n-1

    
        if i >= 2 && i <= m-2
            if Fe(i,j) > 0
                phiU = phi(i-1,j);  
                phiC = phi(i,j);     
                phiD = phi(i+1,j);   

                phi_uds   = phiC;
                phi_smart = smart_face(phiU,phiC,phiD);

                bdc(i,j) = bdc(i,j) - Fe(i,j)*(phi_smart - phi_uds);
            elseif Fe(i,j)<0 && i<=m-3
                phiU = phi(i+2,j);
                phiC = phi(i+1,j);
                phiD = phi(i,j);

                phi_uds   = phiC;
                phi_smart = smart_face(phiU,phiC,phiD);

                bdc(i,j) = bdc(i,j) - Fe(i,j)*(phi_smart - phi_uds);
            end
        end

       
        if j >= 2 && j <= n-2
            if Fn(i,j) > 0
                phiU = phi(i,j-1);   
                phiC = phi(i,j);    
                phiD = phi(i,j+1);   

                phi_uds   = phiC;
                phi_smart = smart_face(phiU,phiC,phiD);

                bdc(i,j) = bdc(i,j) - Fn(i,j)*(phi_smart - phi_uds);
            elseif Fn(i,j) <0 && j<=n-3
                phiU = phi(i,j+2);
                phiC = phi(i,j+1);
                phiD = phi(i,j);

                phi_uds   = phiC;
                phi_smart = smart_face(phiU,phiC,phiD);

                bdc(i,j) = bdc(i,j) - Fn(i,j)*(phi_smart - phi_uds);
            end
        end
        
        if i >= 2 && i <= m-2
            if Fw(i,j) > 0
                phiU = phi(i+1,j);
                phiC = phi(i,j);
                phiD = phi(i-1,j);

                phi_uds   = phiC;
                phi_smart = smart_face(phiU,phiC,phiD);

                bdc(i,j) = bdc(i,j) - Fw(i,j)*(phi_smart - phi_uds);
            elseif Fw(i,j) < 0 && i >=3
                phiU = phi(i-2,j);
                phiC = phi(i-1,j);
                phiD = phi(i,j);

                phi_uds   = phiC;
                phi_smart = smart_face(phiU,phiC,phiD);

                bdc(i,j) = bdc(i,j) - Fw(i,j)*(phi_smart - phi_uds);
            end
        end
        if j>=2 && j<=n-2
            if Fs(i,j)>0
                phiU =phi(i,j+1);
                phiC = phi(i,j);
                phiD = phi(i,j-1);

                phi_uds= phiC;
                phi_smart = smart_face(phiU, phiC, phiD);
                bdc(i,j) = bdc(i,j) - Fs(i,j) * (phi_smart - phi_uds);
            elseif Fs(i,j) < 0 && j >=3
                phiU = phi(i,j-2);
                phiC = phi(i,j-1);
                phiD = phi(i,j);
                
                phi_uds = phiC;
                phi_smart = smart_face(phiU, phiC, phiD);
                bdc(i,j) = bdc(i,j) - Fs(i,j) * (phi_smart - phi_uds);
            end
        end

    end
end
end