%cross_tsp_a.m
%交叉策略A
function cnew=cross_tsp_a(c1,c2,n)
j1=ceil(n*rand)%ceil为向上取整函数
j2=ceil(n*rand)
j3=min(j1,j2)
j4=max(j1,j2)
for j=j3:j4
    i=find(c1==c2(j));
    for k=i:n-1
        c1(k)=c1(k+1);
    end
end
mm=n-(j4-j3+1)+1
j5=ceil(mm*rand)
cnew(1:j5-1)=c1(1:j5-1);
cnew(j5:j5+j4-j3)=c2(j3:j4);
cnew(j5+j4-j3+1:n)=c1(j5:mm-1);
%应该是策略C,在第二个串中随机选择一个交叉区域,加到第一个串的随机的位置。
%删除第一个串中交叉区中出现过的城市。
