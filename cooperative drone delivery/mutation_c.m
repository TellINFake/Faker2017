%mutation_c.m
%±äÒì²ßÂÔC
function cnew=mutation_c(c0,n)
j1=ceil(rand*n);
j2=ceil(rand*n);
j3=min(j1,j2);
j4=max(j1,j2);
cnew=c0;
k=1;
while k<=j4-j3+1
    cnew(j3+k-1)=c0(j4-k+1);
    k=k+1;
end
