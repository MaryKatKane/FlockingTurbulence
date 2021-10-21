function [Kx,Ky,Kthet] = RHS_part(U,V,Om,Xp,Yp,Theta,F)
global P_num;
global P_Phys;

mean_or = nan(P_num.np,1);

Kx = P_Phys.V*cos(Theta) + interp2(P_num.ygrid,P_num.xgrid,[U U(:,1); U(1,:) U(1,1)],Yp,Xp);
Ky = P_Phys.V*sin(Theta) + interp2(P_num.ygrid,P_num.xgrid,[V V(:,1); V(1,:) V(1,1)],Yp,Xp);

Kthet = 0.5*interp2(P_num.ygrid,P_num.xgrid,[Om Om(:,1); Om(1,:) Om(1,1)],Yp,Xp);
Xtot = [Xp; Xp-P_Phys.L; Xp-P_Phys.L; Xp-P_Phys.L; Xp; Xp; Xp+P_Phys.L; Xp+P_Phys.L; Xp+P_Phys.L];
Ytot = [Yp; Yp-P_Phys.L; Yp; Yp+P_Phys.L; Yp-P_Phys.L; Yp+P_Phys.L; Yp-P_Phys.L; Yp; Yp+P_Phys.L];
for ip = 1:P_num.np
    I = mod(find((Xtot-Xp(ip)).^2+(Ytot-Yp(ip)).^2<P_Phys.R^2),P_num.np);
    I(I==0) = P_num.np;
    Kthet(ip) = Kthet(ip)+(mean(mod(Theta(I),2*pi))-Theta(ip))*F;
end

return