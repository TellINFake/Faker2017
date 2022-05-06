% file: ant_pso.m
% 保留每次迭代的最优点
% 以max(t^a*d^(-b))为依据找最优路径，与保留的最优路径比较
 clear all
x=[41 37 54 25 7 2 68 71 54 83 64 18 22 83 91 ...
25 24 58 71 74 87 18 13 82 62 58 45 41 44 4];
y=[94 84 67 62 64 99 58 44 62 69 60 54 60 46 38 ...
38 42 69 71 78 76 40 40 7 32 35 21 26 35 50];
n=30;%n表示城市数目
c=100;%c表示初始信息浓度
q=1000000;
NC=100;
r=0.9;%r表示轨迹持久性
a=1.5;%a表示轨迹相对重要性
b=2;%b表示能见度相对重要性
m=50;%m表示蚂蚁数目
for  i=1:n 
    for j=1:n
dij(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);%dij保存任意两个城市之间的距离
    end
end
for i=1:n
dij(i,i)=0.01;%自己和自己的距离定义为0.01
end
min10=100000;%
t=ones(n)*c;%t表示信息浓度
for i=1:100
    road(i,:)=randperm(n);%产生100条随机路径,randperm(n)表示将1到n随机排列
    ltspc(i)=ca_tsp(n,road(i,:),dij);%记录他们的长度
end
ltspsort=sort(ltspc);%对长度进行排序
l30=ltspsort(m);%选择m个初始解，m只蚂蚁,这里选出第m小的数保存在l30
i1=0;
for i=1:100
    if ltspc(i)<=l30%对这30条较优路径操作
        for k=1:n-1
            t(road(i,k),road(i,k+1))=t(road(i,k),road(i,k+1))+10;%让30条路径留下信息素
            t(road(i,k+1),road(i,k))=t(road(i,k),road(i,k+1));
        end
        t(road(i,1),road(i,n))=t(road(i,1),road(i,n))+10;
        t(road(i,n),road(i,1))=t(road(i,1),road(i,n));
        i1=i1+1;
        pcbest(i1,:)=road(i,:);%初始pcbest,路径,i1是从1-30,等于只重新单独保存这30条路径
        plbest(i1)=ltspc(i);%初始plbest，路径长度
    end
end
[glbest,j]=min(plbest); %初始glbest,路径
gcbest=pcbest(j,:);%初始gcbest,路径
for nc=1:NC
    tabu=ones(m,n);
    tabu(:,1)=0;
    path=ones(m,n);
    for k=1:m
        for step=1:n-1
            ta=t^a;
            tb=dij.^(-b);
            td=ta.*tb;
            pd=tabu(k,:).*td(path(k,step),:);
            pk=pd/sum(pd);
            rk=rand;
            spk=0;
            j=1;
            while j<=n;
                if rk<spk+pk(j)
                    break;
                else
                    spk=spk+pk(j);
                    j=j+1;
                end
            end
            tabu(k,j)=0;
            path(k,step+1)=j;
        end
    end
    for i=1:m
        ltsp(i)=ca_tsp(n,path(i,:),dij);%路径长度
        %交叉操作
        %四种交叉程序分别为 cross_tso_a cross_tso_b cross_tso_c cross_tso_d
        path1(i,:)=cross_tsp_b(path(i,:),gcbest,n);
        path1(i,:)=cross_tsp_b(path1(i,:),pcbest(i,:),n);
        %变异操作
        %四种变异子程序分别为mutation_a mutation_b mutation_c mutation_d
        path1(i,:)=mutation_b(path1(i,:),n);
        %计算新路径长度之和
        ltsp1(i)=ca_tsp(n,path1(i,:),dij);
        %求个体极值
        if ltsp1(i)<ltsp(i)%比较交叉前后的路径长度,选优者保留
            ltsp(i)=ltsp1(i);%ltsp为蚂蚁走过一轮的路径长度,ltsp1为交叉变异后的路径长度
            path(i,:)=path1(i,:);
        end
        if ltsp(i)<plbest(i)%比较选优后的路径
            plbest(i)=ltsp(i);%plbest为初始后的30条路径长度
            pcbest(i,:)=path(i,:);%pcbest为初始后的30条路径
        end
    end
    %找全局极值
    [glbest,j]=min(plbest);
    gcbest=pcbest(j,:);%找到全局最优值
    dt=zeros(n);
    for i=1:m
        for k=1:n-1
            dt(path(i,k),path(i,k+1))=dt(path(i,k),path(i,k+1))+q/ltsp(i);
            dt(path(i,k+1),path(i,k))=dt(path(i,k),path(i,k+1));
        end
        dt(path(i,n),path(i,1))=dt(path(i,n),path(i,1))+q/ltsp(i);
        dt(path(i,1),path(i,n))=dt(path(i,n),path(i,1));
    end
    t=r*t+dt;
end
ta=t.^a;
tb=dij.^(-b);
k=3;
ts(1)=1;
td(:,1)=0;
[my,i]=max(td(1,:));
ts(2)=i;
td(:,i)=0;
while k<=n
    [my,i]=max(td(i,:));
    ts(k)=i;
    td(:,i)=0;
    k=k+1;
end
ltsp0=ca_tsp(n,ts,dij);%ltsp0完全是纯蚁群算法求出来的路径长度
if glbest<ltsp0%glbest是在蚂蚁具有粒子性质的情况下求的的最优路径长度,
     ts=gcbest;%这里的比较相当的重要,通过实验可以知道ltsp0(564.7881)比glbest(426.5438)大很多
     ltsp0=glbest;%蚂蚁的粒子性和初始的随机解很好的保证了全局性,比较好的避免纯蚁群算法容易陷入局部解的缺点
end
k=1;
while k<=n
     x1(k)=x(ts(k));
     y1(k)=y(ts(k));
     k=k+1;
end
x1(n+1)=x1(1);
y1(n+1)=y1(1);
line(x1,y1)
hold on
plot(x,y,'o');
ltsp0
       
        