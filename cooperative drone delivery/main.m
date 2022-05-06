%��ʼ��
clear;
tic;%�����������ʱ��
t1=clock;
alpha=1; %��Ϣ����Ҫ�̶ȵĲ���
beta=5; %����ʽ������Ҫ�̶ȵĲ��� 
rho=0.5; %��Ϣ������ϵ��
max=100; %����������
q=100; %��Ϣ������ǿ��ϵ��
cityNum=50;  %����Ĺ�ģ�����и�����
[dislist,Clist]=tsp(cityNum);
m=cityNum; %���ϸ���
Eta=1./dislist;%EtaΪ�������ӣ�������Ϊ����ĵ���
Tau=ones(cityNum,cityNum);%TauΪ��Ϣ�ؾ���
Tabu=zeros(m,cityNum);%�洢����¼·��������
NC=1;%����������
R_best=zeros(max,cityNum); %�������·��
L_best=inf.*ones(max,1);%�������·�ߵĳ���
L_ave=zeros(max,1);%����·�ߵ�ƽ������

figure(1);
while NC<=max %ֹͣ����֮һ���ﵽ����������
    %��mֻ���Ϸŵ�cityNum��������
    Randpos=[];
    for i=1:(ceil(m/cityNum))
        Randpos=[Randpos,randperm(cityNum)];
    end
    Tabu(:,1)=(Randpos(1,1:m))';
    
    %mֻ���ϰ����ʺ���ѡ����һ�����У���ɸ��Ե�����
    for j=2:cityNum
        for i=1:m
            visited=Tabu(i,1:(j-1)); %�ѷ��ʵĳ���
            J=zeros(1,(cityNum-j+1));%�����ʵĳ���
            P=J;%�����ʳ��е�ѡ����ʷֲ�
            Jc=1;
            for k=1:cityNum
                if length(find(visited==k))==0
                    J(Jc)=k;
                    Jc=Jc+1;
                end
            end
            %�����ѡ���еĸ��ʷֲ�
            for k=1:length(J)
                P(k)=(Tau(visited(end),J(k))^alpha)*(Eta(visited(end),J(k))^beta);
            end
            P=P/(sum(P));
            %������ԭ��ѡȡ��һ������
            Pcum=cumsum(P);
            Select=find(Pcum>=rand);
            to_visit=J(Select(1));
            Tabu(i,j)=to_visit;
        end
    end
    if NC>=2
        Tabu(1,:)=R_best(NC-1,:);
    end
    %��¼���ε������·��
    L=zeros(m,1);
    for i=1:m
        R=Tabu(i,:);
        L(i)=CalDist(dislist,R);
    end
    L_best(NC)=min(L);
    pos=find(L==L_best(NC));
    R_best(NC,:)=Tabu(pos(1),:);
    L_ave(NC)=mean(L);
    drawTSP(Clist,R_best(NC,:),L_best(NC),NC,0);
    NC=NC+1;
    %������Ϣ��
    Delta_Tau=zeros(cityNum,cityNum);
    for i=1:m
        for j=1:(cityNum-1)
            Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+q/L(i);
        end
        Delta_Tau(Tabu(i,cityNum),Tabu(i,1))=Delta_Tau(Tabu(i,cityNum),Tabu(i,1))+q/L(i);
    end
    Tau=(1-rho).*Tau+Delta_Tau;
    Tabu=zeros(m,cityNum); %���ɱ�����
    %pause;
end

%������
Pos=find(L_best==min(L_best));
Shortest_Route=R_best(Pos(1),:);
Shortest_Length=L_best(Pos(1));
figure(2);
plot([L_best L_ave]);
legend('��̾���','ƽ������');
disp(['����������ʱ�䣺',num2str(etime(clock,t1))]);