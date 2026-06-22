%ac*phic+ae*phie+an*phin+as*phis+aw*phiw=bc
%a=-gamma*mag(S)/mag(C)
%aX=b
%element 1: ac*phi1+ae*phi2+an*phi4+0+0=Q1V1(number of phi depends on
%number  of cell containing phi at said location)
%element 2: ac*phi2+ae*phi3+an*phi5+aw*phi1+0=Q2V2
%b=[Q1V1 Q2V2 ... QnVn]
%X=[phi1 phi2 ... phin]
%a has row indexes [phi1 phi2 phi3 phi4 ... phin]
%element1 [ac ae 0 an 0 0 0 0 0]
%element2 [aw ac ae 0 aw 0 0 0 0]
%element3 [0 aw ac 0 0 an 0 0 0]
%...