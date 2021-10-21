function Ek=GetSpectra(Psi,sp1d_d)

global P_num;
vx=1i*P_num.ky.*Psi;
vy=-1i*P_num.kx.*Psi;

%Compute velocity spectrum
Ek=.5*spec1d(vx);
Ek=Ek+.5*spec1d(vy);

end