function [Kthet,U,V]=RHS_omeg(Om)
global P_num;
Psihat = calc_psi(Om);
U = ifft2(1i*P_num.ky.*Psihat,'symmetric');
V = -ifft2(1i*P_num.kx.*Psihat,'symmetric');
Kthet = -fft2(U.*ifft2(1i*P_num.kx.*Om,'symmetric')) ...
    -fft2(V.*ifft2(1i*P_num.ky.*Om,'symmetric'));
return


