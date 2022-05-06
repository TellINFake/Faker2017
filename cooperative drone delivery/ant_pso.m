% file: ant_pso.m
% ����ÿ�ε��������ŵ�
% ��max(t^a*d^(-b))Ϊ����������·�����뱣��������·���Ƚ�
 clear all
x=[41 37 54 25 7 2 68 71 54 83 64 18 22 83 91 ...
25 24 58 71 74 87 18 13 82 62 58 45 41 44 4];
y=[94 84 67 62 64 99 58 44 62 69 60 54 60 46 38 ...
38 42 69 71 78 76 40 40 7 32 35 21 26 35 50];
n=30;%n��ʾ������Ŀ
c=100;%c��ʾ��ʼ��ϢŨ��
q=1000000;
NC=100;
r=0.9;%r��ʾ�켣�־���
a=1.5;%a��ʾ�켣�����Ҫ��
b=2;%b��ʾ�ܼ��������Ҫ��
m=50;%m��ʾ������Ŀ
for  i=1:n 
    for j=1:n
dij(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);%dij����������������֮��ľ���
    end
end
for i=1:n
dij(i,i)=0.01;%�Լ����Լ��ľ��붨��Ϊ0.01
end
min10=100000;%
t=ones(n)*c;%t��ʾ��ϢŨ��
for i=1:100
    road(i,:)=randperm(n);%����100�����·��,randperm(n)��ʾ��1��n�������
    ltspc(i)=ca_tsp(n,road(i,:),dij);%��¼���ǵĳ���
end
ltspsort=sort(ltspc);%�Գ��Ƚ�������
l30=ltspsort(m);%ѡ��m����ʼ�⣬mֻ����,����ѡ����mС����������l30
i1=0;
for i=1:100
    if ltspc(i)<=l30%����30������·������
        for k=1:n-1
            t(road(i,k),road(i,k+1))=t(road(i,k),road(i,k+1))+10;%��30��·��������Ϣ��
            t(road(i,k+1),road(i,k))=t(road(i,k),road(i,k+1));
        end
        t(road(i,1),road(i,n))=t(road(i,1),road(i,n))+10;
        t(road(i,n),road(i,1))=t(road(i,1),road(i,n));
        i1=i1+1;
        pcbest(i1,:)=road(i,:);%��ʼpcbest,·��,i1�Ǵ�1-30,����ֻ���µ���������30��·��
        plbest(i1)=ltspc(i);%��ʼplbest��·������
    end
end
[glbest,j]=min(plbest); %��ʼglbest,·��
gcbest=pcbest(j,:);%��ʼgcbest,·��
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
        ltsp(i)=ca_tsp(n,path(i,:),dij);%·������
        %�������
        %���ֽ������ֱ�Ϊ cross_tso_a cross_tso_b cross_tso_c cross_tso_d
        path1(i,:)=cross_tsp_b(path(i,:),gcbest,n);
        path1(i,:)=cross_tsp_b(path1(i,:),pcbest(i,:),n);
        %�������
        %���ֱ����ӳ���ֱ�Ϊmutation_a mutation_b mutation_c mutation_d
        path1(i,:)=mutation_b(path1(i,:),n);
        %������·������֮��
        ltsp1(i)=ca_tsp(n,path1(i,:),dij);
        %����弫ֵ
        if ltsp1(i)<ltsp(i)%�ȽϽ���ǰ���·������,ѡ���߱���
            ltsp(i)=ltsp1(i);%ltspΪ�����߹�һ�ֵ�·������,ltsp1Ϊ���������·������
            path(i,:)=path1(i,:);
        end
        if ltsp(i)<plbest(i)%�Ƚ�ѡ�ź��·��
            plbest(i)=ltsp(i);%plbestΪ��ʼ���30��·������
            pcbest(i,:)=path(i,:);%pcbestΪ��ʼ���30��·��
        end
    end
    %��ȫ�ּ�ֵ
    [glbest,j]=min(plbest);
    gcbest=pcbest(j,:);%�ҵ�ȫ������ֵ
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
ltsp0=ca_tsp(n,ts,dij);%ltsp0��ȫ�Ǵ���Ⱥ�㷨�������·������
if glbest<ltsp0%glbest�������Ͼ����������ʵ��������ĵ�����·������,
     ts=gcbest;%����ıȽ��൱����Ҫ,ͨ��ʵ�����֪��ltsp0(564.7881)��glbest(426.5438)��ܶ�
     ltsp0=glbest;%���ϵ������Ժͳ�ʼ�������ܺõı�֤��ȫ����,�ȽϺõı��ⴿ��Ⱥ�㷨��������ֲ����ȱ��
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
       
        