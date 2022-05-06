t=[4 5;16 25.8;10 45;20 55;30 65;35 55;29 31;37 26;47 27;30 31.3;31 17;14 7;35.6 13.8;26.7 22.5;21 39;38 42;5 26;28 53;20 13;10 60;26 31;54 38;7 58;12 36;30 2] %24����,��25������origin
save t.mat t
load t.mat
value=[1 1 1 2 3 2 1 3 3 2 2 2 2 2 1 2 3 3 1 1 2 1 1 1]; %24��Ŀ��ļ�ֵ
value=value/100;
time=zeros(1,25); %���UAVʱ�����飬����ŵ��Ƿɻ��ߵĺ��̣������ٶȱ���ʱ�䣬���ٶ�Ϊ��1��
attacktime=zeros(1,25); %���UAVʱ������
%�������������Ӧ���˻��߹��ĺ��̴浽�������һ�����������ʱ�䣬Ȼ�����������ѡ��ĳ����checkһ��ʱ���Ƿ�ϸ񣬺ϸ�Ļ������Դ������������ɱ����ϸ�Ļ���ѡ�θ��ʵ�

%ע��Ŀ���ǰ�����Ŀ��ִ����������������ÿ�ε�������������˻��ջ���ܼ�ֵ��һ������������Ŀ��ļ�ֵ֮�ͣ����Ա�����������ִ�м�ֵ���Ŀ�꣬��ֹ���˻��ɺܾá���ܾú󹥴�Ч�ʱ��
%���������

%%������м��໥����
n=size(t,1);
D=zeros(n,n);
for i=1:n
    for j=1:n
        if i~=j
            D(i,j)=sqrt(sum((t(i,:)-t(j,:)).^2));
        else
            D(i,j)=1e-2;
        end
    end
end
%%��ʼ������
m=10;         %���ϸ���
alpha=1;      %��Ϣ����Ҫ�̶�����
beta=1;       %����������Ҫ�̶�����
gama=2;
rho=0.3;      %��Ϣ�ػӷ�����
Q=1.0;          %����
eta=1./D;     %��������
tau=ones(n,n)+7.1192e-005;%��Ϣ�ؾ���
iter=1;       %����������ʼֵ
iter_max=200;  %�����������ֵ
length_best=zeros(iter_max,1);%ÿ�ε������·������(Ӧ����һ�α�һ��С)
length_ave=zeros(iter_max,1); %ÿ�ε���·��ƽ������ 
%%����Ѱ�����·��
while iter<=iter_max
    whta=cell(8,1);
    lieend=zeros(8,1);
    for zu=1:8
    city_index=1:25;       %��������� 
    table=[];
    start=zeros(4,1);
        temp=randperm(24);
        for i=1:4
        start(i)=temp(i);
        end
    table(:,1)=start;
    j=2;
 while (j<=30)
        for i=1:4
            if i==1 %UAV1ֻ������족����
                if table(1,(j-1))~=25
                    table1=table(1,:);
                    table1=[table1;table(3:4,:)];
                    tabu1=table1(:); %UAV1�Ľ��ɱ������ %25���Ҳ��tabu1��Ļ�����ô
            allow_index1=~ismember(city_index,tabu1);  %���߹��ı��0�����ߵ�Ϊ1������tabu=(1 4)��allow_index=(0 1 1 0 1 1 1...)����ע�⣺allow_index��city_indexͬά��
            allow1=city_index(allow_index1);  %�ѻ����ߵ����ժ�����ˣ������ʵĳ��м��ϣ�
            P1=allow1;
            %������е�ת�Ƹ���
            if numel(allow1)~=0
              for k=1:max(size(allow1))-1
               P1(k)=(tau(table(1,(j-1)),allow1(k))^alpha)*(eta(table(1,(j-1)),allow1(k))^beta)*10000+7.1192e-004;
              end
            P1(max(size(allow1)))=7.1192e-005;
            P1=P1/sum(P1);
            [d1,ind1]=sort(P1,2,'descend');%�Ӵ�С������d1,��Ӧ��ԭ�����ind1
            target1=allow1(ind1(1));
            %���̶ķ�ѡ����һ������
            %pc1=cumsum(P1);  % ��p1 p1+p2 p1+p2+p3 p1+p2+p3+p4 ....����p1<->allow(1)  p2<->allow(2) ...��
            %target_index1=find(pc1>=rand); 
            %target1=allow1(target_index1(1));  %��η��ص���allow�����г��е��������
            table(1,j)=target1;  %��ѡ�������ŵ�·��������
            rr=D(25,table(1,1));
            time(table(1,1))=rr;
            if j>2
            for c=2:(j-1)
                rr=rr+D(table(1,c-1),table(1,c));
            end 
            end
            rrr=rr+D(table(1,j-1),target1);%rrr����UAV1���õ�ʱ�߹��ĺ���
            time(target1)=rrr;
            else
                table(1,j)=25;
            end
                end
                 if table(1,(j-1))==25
                    table(1,j)=25;
                 end
            end          
            if i==2 %UAV2ֻ���𡰴��������
                if (table(2,(j-1))~=25)
                table(2,1)=table(1,1); %�趨����һ�δ������UAV1������Ŀ��
                ta2=table(1:(4*(j-1)+1)); %��ǰԪ��֮ǰ���е�Ԫ��
                tabu21=[];
                tabu22=[];
                tabu2=[];
                for y=1:24
                    if sum(ta2==y)==2
                        tabu21=[tabu21;y];
                    end
                end                       %���ֹ����εķ���tabu21��
                tabu22=setdiff(1:24,ta2); %һ�ζ�û���ֵķ���tabu22��
                tabu2=[tabu21',tabu22];   %tabu2������
            allow_index2=~ismember(city_index,tabu2);  %���߹��ı��0�����ߵ�Ϊ1������tabu=(1 4)��allow_index=(0 1 1 0 1 1 1...)����ע�⣺allow_index��city_indexͬά��
            allow2=city_index(allow_index2);  %�ѻ����ߵ����ժ�����ˣ������ʵĳ��м��ϣ�
            P2=allow2;
            %������е�ת�Ƹ���
           for k=1:(length(allow2)-1)
               P2(k)=tau(table(2,(j-1)),allow2(k))*eta(table(2,(j-1)),allow2(k))*value(allow2(k))*10000;
           end
           P2(max(size(allow2)))=7.1192e-005;
            P2=P2/sum(P2);
            [d2,ind2]=sort(P2,2,'descend');%�Ӵ�С������d1,��Ӧ��ԭ�����ind1
            target2=allow2(ind2(1)); %target2=d1(1);
            %���̶ķ�ѡ����һ������
            %pc2=cumsum(P2);  % ��p1 p1+p2 p1+p2+p3 p1+p2+p3+p4 ....����p1<->allow(1)  p2<->allow(2) ...��
            %target_index2=find(pc2>=rand);  %ѡ���Ǹ����ʽϴ��ѡ�еĵ㣬���ص���allow�����е����
            %target2=allow2(target_index2(1));  %��η��ص���allow�����г��е��������
            %table(2,j)=target2;  %��ѡ�������ŵ�·��������
            oo=D(25,table(2,1));
            attacktime(table(2,1))=oo;
            if j>2
            for c=2:(j-1)
                oo=oo+D(table(2,c-1),table(2,c));
            end 
            end
            ooo=oo+D(table(2,j-1),target2);%ooo����UAV2���õ�ʱ�߹��ĺ���
            if numel(d2)>5
            u=2;
            while (ooo>time(target2)+20 & u<6)
                 target2=allow2(ind2(u));
                 ooo=oo+D(table(2,(j-1)),target2);
                 u=u+1;
            end
            end
            table(2,j)=target2;
            attacktime(target2)=ooo;
                end
                if table(2,(j-1))==25
                    table(2,j)=25;
                end
            end
            if i==3 %UAV3�ǡ��������
                if table(3,(j-1))~=25
                    ta3=table(1:(4*(j-1)+2));
                    tabu3=[];
                    tabu3c=[];
                for y=1:24
                    if sum(ta3==y)==2
                        tabu3=[tabu3;y];
                    end
                end    %�������εķ���tabu3��   
                  for y=1:24
                    if sum(ta3==y)==1
                        tabu3c=[tabu3c;y];
                    end
                end %tabu3c�Ǵ��������������������
            allow_index3=~ismember(city_index,tabu3);  %���߹��ı��0�����ߵ�Ϊ1������tabu=(1 4)��allow_index=(0 1 1 0 1 1 1...)����ע�⣺allow_index��city_indexͬά��
            allow3=city_index(allow_index3);  %�ѻ����ߵ����ժ�����ˣ������ʵĳ��м��ϣ�
            P3=allow3;
            %������е�ת�Ƹ���
           for k=1:(length(allow3)-1)
               %if ismember(allow3(k),tabu3c)==1     
               h=table(3,(j-1))
               P3(k)=(tau(table(3,j-1),allow3(k))^alpha)*(eta(table(3,(j-1)),allow3(k))^beta)*value(allow3(k))*10000+7.1192e-005;%����Ҫ��ģ���Ҫ��ֵ
               %else
               %P3(k)=(tau(table(3,(j-1)),allow3(k))^alpha)*(eta(table(3,(j-1)),allow3(k))^beta)*100+7.1192e-005;%��Щ�Ǵ����ģ�û�м�ֵ
               %end
           end
            P3(max(size(allow3)))=7.1192e-009;
            P3=P3/sum(P3);
            [d3,ind3]=sort(P3,2,'descend');%�Ӵ�С������d1,��Ӧ��ԭ�����ind1
            target3=allow3(ind3(1));
            %���̶ķ�ѡ����һ������
            %pc3=cumsum(P3);  % ��p1 p1+p2 p1+p2+p3 p1+p2+p3+p4 ....����p1<->allow(1)  p2<->allow(2) ...��
            %target_index3=find(pc3>=rand);  %ѡ���Ǹ����ʽϴ��ѡ�еĵ㣬���ص���allow�����е����
            %target3=allow3(target_index3(1));  %��η��ص���allow�����г��е��������
            %table(3,j)=target3;  %��ѡ�������ŵ�·��������
            ww=D(25,table(3,1));
            time(table(3,1))=ww;
            if j>2
            for c=2:(j-1)
                ww=ww+D(table(3,c-1),table(3,c));
            end 
            end
            www=ww+D(table(3,j-1),target3);%www����UAV3���õ�ʱ�߹��ĺ���
            if ismember(target3,tabu3c)==0 %�������
                time(target3)=www;
                table(3,j)=target3;
            else %�������
                attacktime(target3)=www;
                if numel(d3)>5
                u=2;
            while (www>time(target3)+20 & u<6)
                 target3=allow3(ind3(u));
                 www=ww+D(table(3,(j-1)),target3);
                 u=u+1;
            end
                end
            attacktime(target3)=www;
            table(3,j)=target3;%www<time(target3)+10 ˵���˴���������
            end 
                end
                if table(3,(j-1))==25
                    table(3,j)=25;
                end
            end
            if i==4 %UAV4�ǡ��������
                if table(4,(j-1))~=25
                    ta4=table(1:(4*(j-1)+3));
                    tabu4=[];
                    tabu4c=[];
                for y=1:24
                    if sum(ta4==y)==2
                        tabu4=[tabu4;y];
                    end
                end    %�������εķ���tabu4������԰��Ѿ������ķ���tabu4c�У��������ֹ�һ�εģ������ѡ��������tabu4'�е�˵����Ҫ�����Ȼ����һ�����ĺ��̣��ٺ����·���Ƚ� 
                for y=1:24
                    if sum(ta4==y)==1
                        tabu4c=[tabu4c;y];
                    end
                end 
            allow_index4=~ismember(city_index,tabu4);  %���߹��ı��0�����ߵ�Ϊ1������tabu=(1 4)��allow_index=(0 1 1 0 1 1 1...)����ע�⣺allow_index��city_indexͬά��
            allow4=city_index(allow_index4);  %�ѻ����ߵ����ժ�����ˣ������ʵĳ��м��ϣ�
            P4=allow4;
            %������е�ת�Ƹ���
           for k=1:(max(size(allow4))-1)
               %if ismember(allow4(k),tabu4c)==1
               sxx=table(4,(j-1))
               P4(k)=(tau(table(4,(j-1)),allow4(k))^alpha)*(eta(table(4,(j-1)),allow4(k))^beta)*value(allow4(k))*10000+7.1192e-005;
               %else
               %P4(k)=(tau(table(4,(j-1)),allow4(k))^alpha)*(eta(table(4,(j-1)),allow4(k))^beta)*100+7.1192e-005;
               %end
           end
           P4(max(size(allow4)))=7.1192e-009;
            P4=P4/sum(P4);
            [d4,ind4]=sort(P4,2,'descend');%�Ӵ�С������d1,��Ӧ��ԭ�����ind1
            target4=allow4(ind4(1));
            %���̶ķ�ѡ����һ������
            %pc4=cumsum(P4);  % ��p1 p1+p2 p1+p2+p3 p1+p2+p3+p4 ....����p1<->allow(1)  p2<->allow(2) ...��
            %target_index4=find(pc4>=rand);  %ѡ���Ǹ����ʽϴ��ѡ�еĵ㣬���ص���allow�����е����
            %target4=allow4(target_index4(1));  %��η��ص���allow�����г��е��������
            %table(4,j)=target4;  %��ѡ�������ŵ�·��������
            qq=D(25,table(4,1));
            time(table(4,1))=qq;
            if j>2
            for c=2:(j-1)
                qq=qq+D(table(4,c-1),table(4,c));
            end 
            end
            qqq=qq+D(table(4,j-1),target4);%www����UAV3���õ�ʱ�߹��ĺ���
            if ismember(target4,tabu4c)==0 %�������
                time(target4)=qqq;
                table(4,j)=target4;
            else %�������
                attacktime(target4)=qqq;
                if numel(d4)>5
                u=2;
            while (qqq>time(target4)+20 & u<6)
                 target4=allow4(ind4(u));
                 qqq=qq+D(table(4,j-1),target4);
                 u=u+1;
            end
                end
            attacktime(target4)=qqq;
            table(4,j)=target4;%www<time(target4)+10 ˵���˴���������
            end                
                end
                if table(4,(j-1))==25
                    table(4,j)=25;
                end
              end
        end %һ�н���
             if table(:,j)==[25;25;25;25]
                 jishu=0;
                 for i=1:24
                 if sum(table(:)'==i)==2
                     jishu=jishu+1;
                 end
                 end
                 if jishu==24
                     lieend1=j;
                     break;
                 else
                     j=j+1;
                     if j==31
                         break;
                     end
                 end
             else
                 j=j+1;
                 if j==31
                     break;
                 end
             end  
 end %table��4��lieend1�еľ���
 lieend1=j;
 if lieend1==31
     lieend1=lieend1-1;
 end
 lieend(zu)=lieend1;
 disp(['lieend1=' num2str(lieend1)])
 wholetable1=table;
 s1=table;
 disp([wholetable1])
 %��wholetable1�ܹ�4�У���ÿ�н���2-opt֮���ٸ�����Ϣ��
 if lieend1<20
    for tt=1:4
           l1(tt)=0;
                    for e=1:lieend1-1
                        l1(tt)=l1(tt)+D(table(tt,e),table(tt,(e+1)));
                    end
                   % wholetable1(tt,:)=table(tt,:);
        for i=1:lieend1-1
            for j=1:lieend1-1
                if j>i
                    for k=0:fix((j-i)/2)
                        u=s1(tt,(i+k));
                        s1(tt,(i+k))=s1(tt,(j-k));
                        s1(tt,(j-k))=u;
                    end
                    l0=0;
                    for e=1:lieend1-1
                        l0=l0+D(s1(tt,e),s1(tt,(e+1)));
                    end %l0��ÿ��֮��ÿ�еĳ���
                        if (s1(1,1)==25)||(s1(2,1)==25)||(s1(3,1)==25)||(s1(4,1)==25)
                            s1(tt,:)=table(tt,:);
                            break;
                        end
                    if l0<=l1(tt) %%����������������⣬Ӧ���ǰ����ü�ֵ����ָ��ȣ�����ֻ���ú��̴��۱�
                        l1(tt)=l0;
                        wholetable1(tt,:)=s1(tt,:);
                    else
                        l1(tt)=l1(tt);
                        wholetable1(tt,:)=wholetable1(tt,:);
                    end
                    s1(tt,:)=table(tt,:);
                end
            end
        end
    end  %wholetable1��tt,:������ŵĶ��ǸĹ�����̵�·������4�У�l1��tt���ŵ��ǸĹ�����̵ĳ��ȣ���4��
 end
     whta{zu}=wholetable1;
    end %whta�������8��wholetable���ֱ���whta{1} whta{2}��������һ��wholetable��lieend��1���У��ڶ�����lieend��2���С�����
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ܹ�8��С�ֶ�
    %%��ʼ����·����
    
     Length=zeros(8,1);
     for i=1:8
     Length(i)=D(25,whta{i}(1,1))+D(25,whta{i}(2,1))+D(25,whta{i}(3,1))+D(25,whta{i}(4,1));
     end
     for i=1:8
         for j=1:4
             for k=1:(lieend(i)-1)
                 Length(i)=Length(i)+D(whta{i}(j,k),whta{i}(j,k+1));
             end
         end
     end %Length��һ��8�е��ݾ��� ����ÿ��Сwholetable��·������
 J=zeros(8,3);
 for i=1:8
     for j=2:4 %��Ϊ1��������˻�����ֻ����2:4
         for k=1:12 %ֻ����ǰ12������Ϊ���䷽��ÿ������һ��Ҳ����12����������
             if whta{i}(j,k+1)~=whta{i}(j,k);
                 J(i,j-1)=J(i,j-1)+(13-k)*value(whta{i}(j,k))*100; %��ֵ������������ִ�У����Ӧ��ִ�����kС����13-k�ʹ������Ե�Ȩ�ؾʹ󣬵õ�����ָ��J�ʹ�
             end
         end
     end
 end %Jÿһ�е�3�зֱ�װ����ÿ��wholetable�������������У�ÿһ�еļ�ֵ��
 zb=Length'*0.7-sum(J')*0.3; %zhibiao��������Ǽ�ֵ+��������ָ�ꣻ ��8�������о��� %%����length���ָ����ԽСԽ�ã���J'��ֵ���ָ����Խ��Խ�ã���J'��С��length�����������ü���
 
    %%������Сָ��
    if iter==1
        route_best=[];
        [up,upindex]=sort(zb); %up�����Ǵ�С�������кõľ���upindex������������ԭ�����еı�ţ�upindex(1)����̵�wholetable(x)��Ӧ��x
        i=1;
        while(i<=8)
        route_best=whta{upindex(i)};
        if size(route_best,2)==30
            i=i+1;
            route_best=[];
        else
            break;
        end
        end
        length_best(iter)=up(i);
    else
        route_best_last=route_best; %����һ�ε�����route_best��������������һ����
        route_best1=[]; %ÿ�ε������·��  % n:25
        route_best=[];
        [up1,upindex1]=sort(zb);
        i=1;
        while(i<=8)
        route_best1=whta{upindex1(i)};
        if size(route_best1,2)==30
            i=i+1;
            route_best1=[];
        else
            break;
        end
        end
        if up1(i)<=length_best(iter-1)
            route_best=route_best1;
            length_best(iter)=up1(i);
        else
            route_best=route_best_last;
            length_best(iter)=length_best(iter-1);
        end
    end
%%������Ϣ��
%%������Ϣ��
delta_tau=zeros(25,25)+7.1192e-004;
for i=1:8
    for j=1:4
        delta_tau(25,whta{i}(j,1))=Q/zb(i)';%������Ϣ��������ָ����С�ػ��Ǻ�����С��
    end
end
for i=1:8
    for j=1:4
        for k=1:(lieend(i)-1)
            delta_tau(whta{i}(j,k),whta{i}(j,k+1))=delta_tau(whta{i}(j,k),whta{i}(j,k+1))+Q/zb(i)';
        end
    end
end
            
delta_tau=10*delta_tau
tau=(1-rho)*tau+delta_tau
iter=iter+1;
end
%%�����ʾ
[shortest_length,index]=min(length_best);
shortest_route=route_best; %shortest_route ��һ��4�еľ���

disp([ num2str(shortest_route)])
disp(['��̾��룺', num2str(shortest_length)])
shortest_route_1=shortest_route(1,:);
shortest_route_2=shortest_route(2,:);
shortest_route_3=shortest_route(3,:);
shortest_route_4=shortest_route(4,:);

plot([t(25,1);t(shortest_route_1,1);],...
    [t(25,2);t(shortest_route_1,2);],'ro:','LineWidth',1.8);
hold on;
plot([t(25,1);t(shortest_route_2,1);],...
    [t(25,2);t(shortest_route_2,2);],'gp-.','LineWidth',1.8);
hold on;
plot([t(25,1);t(shortest_route_3,1);],...
    [t(25,2);t(shortest_route_3,2);],'cd-','LineWidth',1.8);
hold on;
plot([t(25,1);t(shortest_route_4,1);],...
    [t(25,2);t(shortest_route_4,2);],'m*--','LineWidth',1.8);
hold off;
figure(2)
plot(1:iter_max,length_best,'k')
axis([0,200,100,700])
%legend('��̾���')
xlabel('��������')
ylabel('����ָ��')
%title('������̺�����ƽ�����̱Ƚ�')
%%�������4��UAV�ֱ��ߵĳ���
len1=D(shortest_route_1(1),25);
for s=1:length(shortest_route_1)-1
disp([num2str(len1)])
len1=len1+D(shortest_route_1(s),shortest_route_1(s+1));
end
disp(['UAV 1���̳���', num2str(len1)])
len2=D(shortest_route_2(1),25);
for s=1:length(shortest_route_2)-1
disp([num2str(len2)])
len2=len2+D(shortest_route_2(s),shortest_route_2(s+1));
end
disp(['UAV 2���̳���', num2str(len2)])
len3=D(shortest_route_3(1),25);
for s=1:length(shortest_route_3)-1
disp([num2str(len3)])
len3=len3+D(shortest_route_3(s),shortest_route_3(s+1));
end
disp(['UAV 3���̳���', num2str(len3)])
len4=D(shortest_route_4(1),25);
for s=1:length(shortest_route_4)-1
disp([num2str(len4)])
len4=len4+D(shortest_route_4(s),shortest_route_4(s+1));
end
disp(['UAV 4���̳���', num2str(len4)])
s2=0;
for k=1:12
    s2=s2-value(shortest_route_2(k))*(13-k);
end
s3=0;
for k=1:2:12
    s3=s3-value(shortest_route_3(k))*(13-k);
end
s4=0;
for k=1:2:12
    s4=s4-value(shortest_route_4(k))*(13-k);
end
s=shortest_length/100*0.7+(s2+s3+s4)*0.3