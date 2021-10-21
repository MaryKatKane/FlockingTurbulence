function IIf=getForcingIndex(kinf,ksup,ks,n)
IIf=[];
    for i=1:n
        for j=1:n
            if (sqrt(ks(i,j))>=kinf)&&(sqrt(ks(i,j))<=ksup)
              IIf=[IIf,i+(j-1)*n];
            end
            
        end
    end
end