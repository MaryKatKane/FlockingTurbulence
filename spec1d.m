function Sp=spec1d(field)
global P_num;
Sp=zeros(1,P_num.n/2);
for idk=1:P_num.n/2
    Sp(idk) = sum(real(field(P_num.KK{idk})).^2 ...
        +imag(field(P_num.KK{idk})).^2)/P_num.n^4;
end
return