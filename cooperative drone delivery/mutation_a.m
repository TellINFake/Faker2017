%mutation_a.m
%±äÒì²ßÂÔA
function cnew=mutation_a(c0,n)
j1=ceil(rand*n);
j2=ceil(rand*n);
cnew=c0;
cnew(j1)=c0(j2);
cnew(j2)=c0(j1);
