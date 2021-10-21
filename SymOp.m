function A=SymOp(A)
global P_num
A(1,1)=.5*(A(1,1)+conj(A(1,1)));
A(1,2:P_num.n)=.5*(A(1,2:P_num.n)+conj(fliplr(A(1,2:P_num.n))));
A(2:P_num.n,1)=.5*(A(2:P_num.n,1)+conj(flipud(A(2:P_num.n,1))));
A(2:P_num.n,2:P_num.n)=.5*(A(2:P_num.n,2:P_num.n)+conj(rot90(A(2:P_num.n,2:P_num.n),2)));
end

