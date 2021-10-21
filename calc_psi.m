function Psi = calc_psi(Theta)
global P_num;
Psi = Theta./P_num.ks;
Psi(1) = 0;
return