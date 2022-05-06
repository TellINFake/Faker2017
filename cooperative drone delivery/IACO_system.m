clear
clc

t=[4 5;16 25.8;10 45;20 55;30 65;35 55;29 31;37 26;47 27;50 51.3;51 17;14 7;40 10;26.7 22.5;21 39;38 42;5 26;28 53;20 13;10 60;26 41;54 38;7 58;12 36;30 2]; %24����,��25������origin
save t.mat t
load t.mat

%%������м��໥����
%n=size(t,1);
n=6;
D=zeros(25,25);
for i=1:25
    for j=1:25
        if i~=j
            D(i,j)=sqrt(sum((t(i,:)-t(j,:)).^2));
        else
            D(i,j)=1e-4;
        end
    end
end
%%��ʼ������
m=4;         %���ϸ��� ��ԭ����1ֻ����Ҫ�߹�24���㣬������4ֻ�����ܹ��߹�24���㡿��5�����ϣ�ÿ������4ֻ��ÿֻ�߹�6���㡿
alpha=1;      %��Ϣ����Ҫ�̶�����
beta=5;       %����������Ҫ�̶�����
rho=0.1;      %��Ϣ�ػӷ�����
Q=50;          %����
A=3;
B=1;
eps=1;
lamba=1;
C=3;
D=6;
theta=1;
phai=1.84;
eta=1./D;     %��������
tau=ones(25,25);%��Ϣ�ؾ���
table=zeros(m,n);%·����¼��
iter=1;       %����������ʼֵ
iter_max=400;  %�����������ֵ
route_best=zeros((4*iter_max),n); %ÿ�ε������·��
length_best=zeros(iter_max,1);%ÿ�ε������·������(Ӧ����һ�α�һ��С)
length_ave=zeros(iter_max,1); %ÿ�ε���·��ƽ������

%%����Ѱ�����·��c 
while iter<=iter_max    
    city_index=1:24;       %node ���
    wholetable=[];
    s=1;
    while s<=10
          start=zeros(4,1);
          temp=randperm(24);
            for i=1:4
             start(i)=temp(i);
            end
          table(:,1)=start;
    for j=2:n
        for i=1:m
            tabu=table(1:((j-1)*4+i-1));  %�ѵ�iֻ����֮ǰ���߹�������node������ɱ��� ����table����ӵ�һ������j-1ȫ�ŵ����ɱ��С�
            allow_index=~ismember(city_index,tabu);  %���߹��ı��0�����ߵ�Ϊ1������tabu=(1 4)��allow_index=(0 1 1 0 1 1 1...)����ע�⣺allow_index��city_indexͬά��
            allow=city_index(allow_index);  %�ѻ����ߵ����ժ�����ˣ������ʵĳ��м��ϣ�
            P=allow;
            %������е�ת�Ƹ���
            for k=1:max(size(allow))
                P(k)=(tau(tabu(end-3),allow(k))^alpha)*(eta(tabu(end-3),allow(k))^beta);
                alpha(k)=A*sin(eps*pi*iter./iter_max)+B;
                beta(k)=C*abs(cos(theta*pi*iter./iter_max))+D;
            end 
            P=P/sum(P);
            %���̶ķ�ѡ����һ��node
            pc=cumsum(P);  % ��p1 p1+p2 p1+p2+p3 p1+p2+p3+p4 ....����p1<->allow(1)  p2<->allow(2) ...��
            target_index=Q./length_best;
            tabu(i+1)=tabu(i)+delta_tau+e*delta;
            target_index=find(pc>=rand);  %ѡ���Ǹ����ʽϴ��ѡ�еĵ㣬���ص���allow�����е����
            target=allow(target_index(1));  %��η��ص���allow�����г��е��������
            table(i,j)=target;  %��ѡ�������ŵ�·��������
        end
    end
    wholetable=[wholetable;table]; 
    table=zeros(m,n); %����Ҫ�����ǰ�table������Ȼ��������
    s=s+1;                %��5�������Ժ�table���һ�������������´�iterѭ����ʱ��Ӧ�����㣬��ʵ������Ҳ�У��������ٸ�table��һ�и�����
    end
    %��wholetable�ܹ�40�У���ÿ�н���l1  and t1֮���ٸ�����Ϣ��
     s1=wholetable;
     wholetable2=wholetable;
    for tt=1:40
       % s1(tt,:)=wholetable(tt,:);
           l1(tt)=0;
                    for e=1:5
                        l1(tt)=l1(tt)+D(wholetable(tt,e),wholetable(tt,(e+1)));
                    end
                    l1(tt)=l1(tt)+D(wholetable(tt,1),25)+D(wholetable(tt,6),25); %l1��2opt֮ǰÿ�еĳ���
                   % wholetable2(tt,:)=wholetable(tt,:);
        for i=1:5
            for j=1:6
                if j>i
                    for k=0:fix((j-i)/2)
                        u=s1(tt,(i+k));
                        s1(tt,(i+k))=s1(tt,(j-k));
                        s1(tt,(j-k))=u;
                    end
                    l0=0;
                    for e=1:5
                        l0=l0+D(s1(tt,e),s1(tt,(e+1)));
                    end
                    l0=l0+D(s1(tt,1),25)+D(s1(tt,6),25); %l0��ÿ��֮��ÿ�еĳ���
                    if l0<=l1(tt)
                        l1(tt)=l0;
                        wholetable2(tt,:)=s1(tt,:);
                    else
                        l1(tt)=l1(tt);
                        wholetable2(tt,:)=wholetable2(tt,:);
                    end
                    s1(tt,:)=wholetable(tt,:);
                end
            end
        end
    end  %wholetable2��tt,:������ŵĶ��ǸĹ�����̵�·������40�У�l1��tt���ŵ��ǸĹ�����̵ĳ��ȣ���40��
    %%����·��
    Length_=[];
    h=0;
    while h<=37
    Length1=zeros(4,1);
    for i=1:4
        route1=wholetable2((i+h),:);  %���ѵ�iֻ�����߹���·����route����route�����mֻ����ÿֻ�߹������С�
        for j=1:5
            Length1(i)=Length1(i)+D(route1(j),route1(j+1));
        end
        Length1(i)=Length1(i)+D(route1(6),25)+D(route1(1),25);
    end
    length1=0;
    for r=1:4
        length1=length1+Length1(r);
    end
     Length_=[Length_;length1];  %Length_��һ��10�е��ݾ���
     h=h+4;
    end
    %%�������·����ƽ������
    if iter==1
        [min_length,min_index]=min(Length_);
        length_best(iter)=min_length;
        length_ave(iter)=mean(Length_);
        route_best(1:4,:)=wholetable2((4*min_index-3):(4*min_index),:); %route_best(1:4,:)�Ǳ���10ֻ�������ܺ���̵�4��·��
    else
        [min_length,min_index]=min(Length_);
        length_best(iter)=min(min_length,length_best(iter-1));
        length_ave(iter)=mean(Length_);
        c(i,j)=Ti./(Tj-Ti);
        ita(i,j)=c(i,j)*w./length_best;
        if length_best(iter)==min_length
          route_best((4*iter-3):(4*iter),:)=wholetable2((4*min_index-3):(4*min_index),:);
        else
            route_best((4*iter-3):(4*iter),:)=route_best((4*iter-7):(4*iter-4),:);
        end
    end  
  %%������Ϣ��
delta_tau=zeros(25,25);
for i=1:40
    for j=1:5
        delta_tau(wholetable2(i,j),wholetable2(i,j+1))=delta_tau(wholetable2(i,j),wholetable2(i,j+1))+Q/l1(i);
    end
    delta_tau(wholetable2(i,6),25)=delta_tau(wholetable2(i,6),25)+Q/l1(i);
    delta_tau(25,wholetable2(i,6))=delta_tau(25,wholetable2(i,1))+Q/l1(i);
end
for i=1:k
    rho(k)=phai*Gaussion()./(1+lamda*e^(iter/iter_max));
    tau=(1-rho(k))*tau+delta_tau;
    iter=iter+1;
    table=zeros(m,n);
end
end
    %%�����ʾ
[shortest_length,index]=min(length_best);
shortest_route=route_best((4*index-3):(4*index),:);  %shortest_route ��һ��4�еľ���
shortest_route_1=shortest_route(1,:);
shortest_route_2=shortest_route(2,:);
shortest_route_3=shortest_route(3,:);
shortest_route_4=shortest_route(4,:);
disp(['��̾��룺', num2str(shortest_length)])
disp(['���·����', num2str(shortest_route_1),'      ',num2str(shortest_route_2),'     ',num2str(shortest_route_3),'      ',num2str(shortest_route_4)])
%%�������4��UAV�ֱ��ߵĳ���
len1=0;
for s=1:5
len1=len1+D(shortest_route_1(s),shortest_route_1(s+1));
end
len1=len1+D(shortest_route_1(1),25)+D(shortest_route_1(6),25);
disp(['UAV 1���̳���', num2str(len1)])
len2=0;
for s=1:5
len2=len2+D(shortest_route_2(s),shortest_route_2(s+1));
end
len2=len2+D(shortest_route_2(1),25)+D(shortest_route_2(6),25);
disp(['UAV 2���̳���', num2str(len2)])
len3=0;
for s=1:5
len3=len3+D(shortest_route_3(s),shortest_route_3(s+1));
end
len3=len3+D(shortest_route_3(1),25)+D(shortest_route_3(6),25);
disp(['UAV 3���̳���', num2str(len3)])
len4=0;
for s=1:5
len4=len4+D(shortest_route_4(s),shortest_route_4(s+1));
end
len4=len4+D(shortest_route_4(1),25)+D(shortest_route_4(6),25);
disp(['UAV 4���̳���', num2str(len4)])

%%��ͼ
figure(1)
%plot([t(25,1);t(shortest_route_1,1);t(25,1);t(shortest_route_2,1);t(25,1);t(shortest_route_3,1);t(25,1);t(shortest_route_4,1);t(25,1)],...
%    [t(25,2);t(shortest_route_1,2);t(25,2);t(shortest_route_2,2);t(25,2);t(shortest_route_3,2);t(25,2);t(shortest_route_4,2);t(25,2)],'o-');
plot([t(25,1);t(shortest_route_1(1),1);t(shortest_route_1(2),1);t(shortest_route_1(3),1);t(shortest_route_1(4),1);t(shortest_route_1(5),1);t(shortest_route_1(6),1);t(25,1)],...
    [t(25,2);t(shortest_route_1(1),2);t(shortest_route_1(2),2);t(shortest_route_1(3),2);t(shortest_route_1(4),2);t(shortest_route_1(5),2);t(shortest_route_1(6),2);t(25,2)],'ro:','LineWidth',1.8);
hold on;
plot([t(25,1);t(shortest_route_2(1),1);t(shortest_route_2(2),1);t(shortest_route_2(3),1);t(shortest_route_2(4),1);t(shortest_route_2(5),1);t(shortest_route_2(6),1);t(25,1)],...
    [t(25,2);t(shortest_route_2(1),2);t(shortest_route_2(2),2);t(shortest_route_2(3),2);t(shortest_route_2(4),2);t(shortest_route_2(5),2);t(shortest_route_2(6),2);t(25,2)],'bp-.','LineWidth',1.8);
hold on;
plot([t(25,1);t(shortest_route_3(1),1);t(shortest_route_3(2),1);t(shortest_route_3(3),1);t(shortest_route_3(4),1);t(shortest_route_3(5),1);t(shortest_route_3(6),1);t(25,1)],...
    [t(25,2);t(shortest_route_3(1),2);t(shortest_route_3(2),2);t(shortest_route_3(3),2);t(shortest_route_3(4),2);t(shortest_route_3(5),2);t(shortest_route_3(6),2);t(25,2)],'kd-','LineWidth',1.8);
hold on;
plot([t(25,1);t(shortest_route_4(1),1);t(shortest_route_4(2),1);t(shortest_route_4(3),1);t(shortest_route_4(4),1);t(shortest_route_4(5),1);t(shortest_route_4(6),1);t(25,1)],...
    [t(25,2);t(shortest_route_4(1),2);t(shortest_route_4(2),2);t(shortest_route_4(3),2);t(shortest_route_4(4),2);t(shortest_route_4(5),2);t(shortest_route_4(6),2);t(25,2)],'m*--','LineWidth',1.8);
hold off;
for i=1:24
    text(t(i,1),t(i,2),['  T' num2str(i)]);
end
text(50,150,['\leftarrowUAV 1'])
text(10,550,['UAV 2\rightarrow'])
text(350,180,['\leftarrowUAV 3'])
text(360,60,['\leftarrowUAV 4'])
grid on
%text(t(shortest_route(1),1),t(shortest_route(2),2),'   ���');
%text(t(shortest_route(end),1),t(shortest_route(end),2),'    �յ�');
xlabel('X/m')
ylabel('Y/m')
title(['Cooperative task scheduling result of IACO'])

figure(2)
plot(t(:,1),t(:,2),'k*');
text(300,20,['delivery center'])
text(360,480,['terminal customer'])
xlabel('X/m');
ylabel('Y/m');
hold off;
