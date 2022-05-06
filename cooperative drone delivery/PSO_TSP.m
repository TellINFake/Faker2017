function varargout = PSO_TSP(varargin)

t=[40 50;160 258;100 450;200 550;300 650;350 550;290 310;370 260;470 270;500 513;510 170;140 70;400 100;267 225;210 390;380 420;50 260;280 530;200 130;100 600;260 410;540 380;70 580;120 360;300 20]; %24个点,第25个点事origin
save t.mat t
load t.mat

%初始化
Alpha=0.25; %个体经验保留概率
Beta=0.25; %全局经验保留概率 
NC_max=500; %最大迭代次数
m=20;  %微粒数
CustomNum=24;  %问题的规模（客户个数）
[dislist,Clist]=tsp(CustomNum);

NC=1;%迭代计数器
R_best=zeros(NC_max,CustomNum); %各代最佳路线
L_best=inf.*ones(NC_max,1);%各代最佳路线的长度
L_ave=zeros(NC_max,1);%各代路线的平均长度

%产生微粒的初始位置
for i=1:m
    x(i,:)=randperm(CustomNum);
    L(i)=CalDist(dislist,x(i,:));
end
p=x; %p为个体最好解
pL=L;
[L_best(1), n_best]=min(L);
R_best(1,:)=x(n_best,:);
L_ave(1,1)=mean(L);

%初始交换序
v=ones(CustomNum-1,2,m)*(round(rand*(CustomNum-1))+1);

figure(1);
while NC<=NC_max %停止条件之一：达到最大迭代次数
    for i=1:m
        xnew(i,:)=changeFun(x(i,:),v(:,:,i));
        A=changeNum(x(i,:),p(i,:));
        Arand=randFun(A,Alpha);
        xnew(i,:)=changeFun(xnew(i,:),Arand);
        B=changeNum(x(i,:),R_best(NC,:));
        Brand=randFun(B,Beta);
        xnew(i,:)=changeFun(xnew(i,:),Brand);
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
    drawTSP10(Clist,R_best(NC,:),L_best(NC,1),NC,0);
    %pause;
    NC=NC+1;
end

%输出结果
Pos=find(L_best==min(L_best));
Shortest_Route=R_best(Pos(1),:);
Shortest_Length=L_best(Pos(1));
figure(2);
plot([L_best L_ave]);
legend('最短距离','平均距离');
end

function xnew=changeFun(x,C);
changeLen=size(C,1);
xnew=x;
for i=1:changeLen
    a=xnew(C(i,1));
    xnew(C(i,1))=xnew(C(i,2));
    xnew(C(i,2))=a;
end
end

function C=changeNum(x,y);
CustomNum=size(x,2);
C=ones(CustomNum-1,2);
for i=1:CustomNum-1
    pos=find(x==y(i));
    C(i,:)=[i pos];
    x=changeFun(x,C(i,:));
end
end

function v=randFun(v,w);
randLen=size(v,1);

for i=1:randLen
    if rand>w
        v(i,2)=v(i,1);
    end
end    
end