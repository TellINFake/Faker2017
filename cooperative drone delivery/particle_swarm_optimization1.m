clear

%初始化
Alpha=0.85; %个体经验保留概率
Beta=0.85; %全局经验保留概率 
NC_max=2000; %最大迭代次数
m=30;  %微粒数
CityNum=30;  %问题的规模（城市个数）
[dislist,Clist]=tsp(CityNum);

NC=1;%迭代计数器
R_best=zeros(NC_max,CityNum); %各代最佳路线
L_best=inf.*ones(NC_max,1);%各代最佳路线的长度
L_ave=zeros(NC_max,1);%各代路线的平均长度

%产生微粒的初始位置
for i=1:m
    x(i,:)=randperm(CityNum);
    L(i)=CalDist(dislist,x(i,:));
end
p=x; %p为个体最好解
pL=L;
[L_best(1,1) n_best]=min(L);
R_best(1,:)=x(n_best,:);
L_ave(1,1)=mean(L);

%初始交换序
v=ones(CityNum-1,2,m);

figure(1);
while NC<=NC_max %停止条件之一：达到最大迭代次数
    for i=1:m
        xnew(i,:)=changeFun(x(i,:),v(:,:,i));
        if rand<Alpha
            A=changeNum(x(i,:),p(i,:));
            xnew(i,:)=changeFun(xnew(i,:),A);
        end
        if rand<Beta
            B=changeNum(x(i,:),R_best(NC,:));
            xnew(i,:)=changeFun(xnew(i,:),B);
        end
        v(:,:,i)=changeNum(x(i,:),xnew(i,:));
        L(i)=CalDist(dislist,xnew(i,:));
        if L(i)<pL(i)
            p(i,:)=xnew(i,:);
            pL(i)=L(i);
        end
    end
    [L_bestnew n_best]=min(L);
    R_bestnew=xnew(n_best,:);
    L_ave(NC+1,1)=mean(L);
    if L_bestnew<L_best(NC,1)
        L_best(NC+1,1)=L_bestnew;
        R_best(NC+1,:)=R_bestnew;
    else
        L_best(NC+1,1)=L_best(NC,1);
        R_best(NC+1,:)=R_best(NC,:);
    end
    x=xnew;
    drawTSP(Clist,R_best(NC,:),L_best(NC,1),NC,0);
    NC=NC+1;
end

%输出结果
Pos=find(L_best==min(L_best));
Shortest_Route=R_best(Pos(1),:);
Shortest_Length=L_best(Pos(1));
figure(2);
plot([L_best L_ave]);
legend('最短距离','平均距离');
