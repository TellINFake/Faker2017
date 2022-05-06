%mutation_d.m
%±äÒì²ßÂÔD
function cnew=mutation_d(c0,n)
j1=ceil(rand*n);
j2=ceil(rand*n);
j3=min(j1,j2);
j4=max(j1,j2);
cnew=c0;
if j4>j3
    k=1;
    while k<=j4-j3
        cnew(j3+k-1)=c0(j3+k);
        k=k+1;
    end
    cnew(j4)=c0(j3);
end
