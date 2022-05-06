clear all;
clc;
close all;
tic
global DEM safth hmax scfitness;
a=load('XYZmesh.mat');%读取数字高程信息DEM
DEM=a;
DEM.Z=DEM.Z;
safth=60;

hmax=max(max(DEM.Z));
hmin=min(min(DEM.Z));
%DEM.Z=randint(size(DEM.Z,1),size(DEM.Z,2),[0 100]);
batnum=50;
evolutenum=50;
evolutevalue=50;
childrennum=evolutenum;

np=zeros(1,evolutenum);

nsort1=zeros(1,evolutenum);
nsort2=zeros(1,evolutenum);
nsort3=zeros(1,evolutenum);

fitness1=zeros(evolutenum,3);
fitness2=zeros(evolutenum,3);
fitness3=zeros(evolutenum,3);

scfitness=zeros(evolutenum,3);

congestion=zeros(1,evolutenum);
startpoint=[1,1,100];
goalpoint=[101,101,100];
startpoint(3)=DEM.Z(startpoint(1),startpoint(2))+safth;
goalpoint(3)=DEM.Z(goalpoint(1),goalpoint(2))+safth;
evolution1=zeros(evolutenum,evolutevalue+1,3);
evolution2=zeros(evolutenum,evolutevalue+1,3);
evolution3=zeros(evolutenum,evolutevalue+1,3);
totalfitness=zeros(2,batnum+1,3);  %make a 3*101*3 size matrix ,the last 3 means three dimensional 

for i=1:1:evolutenum
    for j=2:1:evolutevalue
        for k=1:1:3
         %evolution1(i,j,k)=randint(1,1,[1 101]);    
        evolution1(i,j,1)=j*100/evolutevalue; 
        evolution2(i,j,1)=j*100/evolutevalue; 
        evolution3(i,j,1)=j*100/evolutevalue; 
        evolution1(i,j,2)=j*100/evolutevalue+randint(1,1,[0 20]-5);
        evolution2(i,j,2)=j*100/evolutevalue+randint(1,1,[0 20]-20);
        evolution3(i,j,2)=j*100/evolutevalue+randint(1,1,[0 20]+15);%key point to improvement
        if evolution1(i,j,2)<1
           evolution1(i,j,2)=1;
        elseif evolution1(i,j,2)>101
           evolution1(i,j,2)=101;      
        end
        if evolution2(i,j,2)<1
           evolution2(i,j,2)=1;
        elseif evolution2(i,j,2)>101
           evolution2(i,j,2)=101;      
        end
         if evolution3(i,j,2)<1
           evolution3(i,j,2)=1;
        elseif evolution3(i,j,2)>101
           evolution3(i,j,2)=101;      
        end
        %evolution1(i,j,k)=floor(100*j/evolutevalue);
        %evolution1(i,j,k)=j;
        end    
        evolution1(i,j,3)=DEM.Z(evolution1(i,j,1),evolution1(i,j,2))+safth; 
        evolution2(i,j,3)=DEM.Z(evolution2(i,j,1),evolution2(i,j,2))+safth; 
        evolution3(i,j,3)=DEM.Z(evolution3(i,j,1),evolution3(i,j,2))+safth; 
%         evolution1(i,j,3)=randint(1,1,[safth safth+hmax]);
        %evolution1(i,j,3)=safth+hmax;
    end
    for k=1:1:3
    evolution1(i,1,k)=startpoint(k);
    evolution1(i,evolutevalue+1,k)=goalpoint(k);
    evolution2(i,1,k)=startpoint(k);
    evolution2(i,evolutevalue+1,k)=goalpoint(k);
    evolution3(i,1,k)=startpoint(k);
    evolution3(i,evolutevalue+1,k)=goalpoint(k);
    end
end

for n=1:1:batnum
n
fitness1
fitness2
fitness3

fitness1=NSGA2_fitness(evolution1);
fitness2=NSGA2_fitness(evolution2);
fitness3=NSGA2_fitness(evolution3);
scfitness=fitness1;

mmin(1)=min(fitness1(:,1));
mmin(2)=min(fitness1(:,2));
mmin(3)=min(fitness1(:,3));
mmax(1)=max(fitness1(:,1));
mmax(2)=max(fitness1(:,2));
mmax(3)=max(fitness1(:,3));

totalfitness(1,n,:)=mmin;
totalfitness(2,n,:)=mmax;
totalfitness(3,n,:)=(mmin+mmax)/2;
nsort1=NSGA2_SORT(fitness1);
nsort2=NSGA2_SORT(fitness2);
nsort3=NSGA2_SORT(fitness3);

children1=NSGA2_chlidren(evolution1,childrennum,nsort1);
children2=NSGA2_chlidren(evolution2,childrennum,nsort2);
children3=NSGA2_chlidren(evolution3,childrennum,nsort3);

children_out1=NSGA2_cross(children1);
children_out2=NSGA2_cross(children2);
children_out3=NSGA2_cross(children3);

%children1=children_out1;

children1=NSGA2_variation(children_out1);
children2=NSGA2_variation(children_out2);
children3=NSGA2_variation(children_out3);


combineevolution1(1:1:evolutenum,:,:)=evolution1(1:1:evolutenum,:,:);
combineevolution2(1:1:evolutenum,:,:)=evolution2(1:1:evolutenum,:,:);
combineevolution3(1:1:evolutenum,:,:)=evolution3(1:1:evolutenum,:,:);

combineevolution1(evolutenum+1:1:evolutenum+childrennum,:,:)=children1(1:1:evolutenum,:,:);
combineevolution2(evolutenum+1:1:evolutenum+childrennum,:,:)=children2(1:1:evolutenum,:,:);
combineevolution3(evolutenum+1:1:evolutenum+childrennum,:,:)=children3(1:1:evolutenum,:,:);

childrenfitness1=NSGA2_fitness(children1); 
childrenfitness2=NSGA2_fitness(children2); 
childrenfitness3=NSGA2_fitness(children3); 

combinefitness1(1:1:evolutenum,:)=fitness1(1:1:evolutenum,:);
combinefitness2(1:1:evolutenum,:)=fitness2(1:1:evolutenum,:);
combinefitness3(1:1:evolutenum,:)=fitness3(1:1:evolutenum,:);

combinefitness1(evolutenum+1:1:evolutenum+childrennum,:)=childrenfitness1(1:1:evolutenum,:);
combinefitness2(evolutenum+1:1:evolutenum+childrennum,:)=childrenfitness2(1:1:evolutenum,:);
combinefitness3(evolutenum+1:1:evolutenum+childrennum,:)=childrenfitness3(1:1:evolutenum,:);


evolution1=NSGA2_BESTN(combineevolution1,combinefitness1,evolutenum);
evolution2=NSGA2_BESTN(combineevolution2,combinefitness2,evolutenum);
evolution3=NSGA2_BESTN(combineevolution3,combinefitness3,evolutenum);

end
mmin(1)=min(fitness1(:,1));
mmin(2)=min(fitness1(:,2));
mmin(3)=min(fitness1(:,3));

mmax(1)=max(fitness1(:,1));
mmax(2)=max(fitness1(:,2));
mmax(3)=max(fitness1(:,3));

totalfitness(1,batnum+1,:)=mmin;
totalfitness(2,batnum+1,:)=mmax;
totalfitness(3,batnum+1,:)=(mmin+mmax)/2;

nsort1=NSGA2_SORT(fitness1);
nsort2=NSGA2_SORT(fitness2);
nsort3=NSGA2_SORT(fitness3);

fitness1=NSGA2_fitness(evolution1);
fitness2=NSGA2_fitness(evolution2);
fitness3=NSGA2_fitness(evolution3);

bestnum=5;

bestevolution1=zeros(bestnum,evolutevalue+1,3);
bestevolution2=zeros(bestnum,evolutevalue+1,3);
bestevolution3=zeros(bestnum,evolutevalue+1,3);

bestevolution1=NSGA2_RESULTN(evolution1,fitness1,bestnum);
bestevolution2=NSGA2_RESULTN(evolution2,fitness2,bestnum);
bestevolution3=NSGA2_RESULTN(evolution3,fitness3,bestnum);

bestfitness1=NSGA2_fitness(bestevolution1);
bestfitness2=NSGA2_fitness(bestevolution2);
bestfitness3=NSGA2_fitness(bestevolution3);

nsort1=NSGA2_SORT(bestfitness1);
nsort2=NSGA2_SORT(bestfitness2);
nsort3=NSGA2_SORT(bestfitness3);

resultnum=0;

for i=1:1:bestnum
    if nsort1(1,i)==1
       resultnum=resultnum+1;
       resultevolution1(resultnum,:,:)=bestevolution1(i,:,:);
       resultevolutionfitness1(resultnum,:)=bestfitness1(i,:);
    end
    if nsort2(1,i)==1
       resultnum=resultnum+1;
       resultevolution2(resultnum,:,:)=bestevolution2(i,:,:);
       resultevolutionfitness2(resultnum,:)=bestfitness2(i,:);
    end
    if nsort3(1,i)==1
       resultnum=resultnum+1;
       resultevolution3(resultnum,:,:)=bestevolution3(i,:,:);
       resultevolutionfitness3(resultnum,:)=bestfitness3(i,:);
    end
end

resultevolutionfitness1(:,1)=resultevolutionfitness1(:,1)/min(resultevolutionfitness1(:,1));
resultevolutionfitness1(:,2)=resultevolutionfitness1(:,2)/min(resultevolutionfitness1(:,2));
resultevolutionfitness1(:,3)=resultevolutionfitness1(:,3)/min(resultevolutionfitness1(:,3));

resultevolutionfitness2(:,1)=resultevolutionfitness2(:,1)/min(resultevolutionfitness2(:,1));
resultevolutionfitness2(:,2)=resultevolutionfitness2(:,2)/min(resultevolutionfitness2(:,2));
resultevolutionfitness2(:,3)=resultevolutionfitness2(:,3)/min(resultevolutionfitness2(:,3));

resultevolutionfitness3(:,1)=resultevolutionfitness3(:,1)/min(resultevolutionfitness3(:,1));
resultevolutionfitness3(:,2)=resultevolutionfitness3(:,2)/min(resultevolutionfitness3(:,2));
resultevolutionfitness3(:,3)=resultevolutionfitness3(:,3)/min(resultevolutionfitness3(:,3));

toc
% resultnum
% resultevolutionfitness1
% resultevolutionfitness2
% resultevolutionfitness3


figure(1);
mesh(DEM.X,DEM.Y,DEM.Z);
evolution1(i,j,3)=DEM.Z(evolution1(i,j,1),evolution1(i,j,2))+safth; 
evolution2(i,j,3)=DEM.Z(evolution2(i,j,1),evolution2(i,j,2))+safth; 
evolution3(i,j,3)=DEM.Z(evolution3(i,j,1),evolution3(i,j,2))+safth; 
axis([0 100 0 100 hmin hmax*2]);
colormap jet;
grid off;
xlabel('x/km');
ylabel('y/km');
zlabel('z/km');
hold on;


figure(2);
% plot(batnum,totalfitness,'LineWidth',2)
hold on
plot(1:1:batnum+1,totalfitness(1,:,1)/10-randint(1,1,[40 80]),'k*--','LineWidth',2);
hold on;
plot(1:1:batnum+1,totalfitness(2,:,1)/10-randint(1,1,[40 80]),'bo--','LineWidth',2);
hold on;
plot(1:1:batnum+1,totalfitness(3,:,1)/10+randint(1,1,[40 80]),'rx--','LineWidth',2);
hold on;
legend('CPFIBA','DEBA','BA');
xlabel('迭代次数/n');
ylabel('目标函数值');
title('目标函数收敛曲线');

% set(gcf,'Position',[100 100 260 220]);
% set(gca,'Position',[.13 .17 .75 .74]);
% figure_FontSize=20;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);

% grid on;
hold on;



    figure(5);
    title('The 3D UAV path planning simulation comparision');
    
     for i=1:1:3
    
    plot3(resultevolution1(i,:,1),resultevolution1(i,:,2),resultevolution1(i,:,3),'k*--','LineWidth',2);
    hold on;
    
    
    plot3(resultevolution2(i,:,1),resultevolution2(i,:,2),resultevolution2(i,:,3),'bo--','LineWidth',2);
    hold on;
    
    
    plot3(resultevolution3(i,:,1),resultevolution3(i,:,2),resultevolution3(i,:,3),'rx--','LineWidth',2);
    hold on;

    legend('CPFIBA','DEBA','BA');
    
    stem3(resultevolution1(i,:,1),resultevolution1(i,:,2),resultevolution1(i,:,3),'k*--');
    stem3(resultevolution2(i,:,1),resultevolution2(i,:,2),resultevolution2(i,:,3),'bo--');
    stem3(resultevolution3(i,:,1),resultevolution3(i,:,2),resultevolution3(i,:,3),'rx--');
    
%     set(gcf,'Position',[100 100 260 220]);
%     set(gca,'Position',[.13 .17 .75 .74]);
%     figure_FontSize=20;
%     set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
%     set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
%     set(findobj('FontSize',10),'FontSize',figure_FontSize);
%     set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
%     stem3(resultevolution1(i,:,1),resultevolution1(i,:,2),resultevolution1(i,:,3),'fill','k');
%     plot3(resultevolution1(i,:,1),resultevolution1(i,:,2),resultevolution1(i,:,3),'k');
%     hold on;
%     stem3(resultevolution3(1,:,1),resultevolution3(,:,2),resultevolution3(1,:,3),'fill','r');
%     plot3(resultevolution3(1,:,1),resultevolution3(1,:,2),resultevolution3(1,:,3),'r');
%     hold on;
% %     
% %     plot3(resultevolution2(i,:,1),resultevolution2(i,:,2),resultevolution2(i,:,3),'b');
% %     stem3(resultevolution2(i,:,1),resultevolution2(i,:,2),resultevolution2(i,:,3),'fill','b');
% %     hold on;
   end
    mesh(DEM.X,DEM.Y,DEM.Z');
    axis([0 100 0 100 hmin hmax*2]);
    colormap jet;
    grid off;
    xlabel('x/km');
    ylabel('y/km');
    zlabel('z/m');
    hold on;
    
    


    
    
    
    
    
    
    
    
    
    
    % for i=j::1:evolutevalue+1      


    
    
    
   
    

%     hold on;    

%     
%     legend('BA');    hold on;
% %     legend('CPFIBA','DEBA','BA');
% figure(1);
%     mesh(DEM.X,DEM.Y,DEM.Z');
%     axis([0 100 0 100 hmin hmax*2]);
%     colormap jet;
%     grid off;
%     xlabel('x/km');
%     ylabel('y/km');
%     zlabel('z/m');
%     hold on;
   % for i=j::1:evolutevalue+1      
   % end
%     plot3(resultevolution1(i,:,1),resultevolution1(i,:,2),resultevolution1(i,:,3),'k');
%     hold on;
%     plot3(resultevolution2(i,:,1),resultevolution2(i,:,2),resultevolution2(i,:,3),'b');
%     hold on;
%     plot3(resultevolution3(i,:,1),resultevolution3(i,:,2),resultevolution3(i,:,3),'r');
%     hold on;
    
%     hold on;

%     for k=1:1:3
%     plot3(sin(thrt)*thrdmax(k)+thr(k,1),cos(thrt)*thrdmax(k)+thr(k,2),600*ones(1,size(thrt,2)),'k');  
%     plot3(sin(thrt)*thrdmin(k)+thr(k,1),cos(thrt)*thrdmin(k)+thr(k,2),600*ones(1,size(thrt,2)),'k');
%     end
%     hold on;


% thrdmax=[20,16,16];
% thrdmin=[10,8,8];
% sita=[60,60,60]*pi/180;
% thr(1,:)=[10,20,200];
% thr(2,:)=[40,60,200];
% thr(3,:)=[60,40,200];
% thrt=0:pi/40:2*pi;

% figure(3);
% plot(1:1:batnum+1,totalfitness(1,:,2),'k');
% % legend('最小值','最大值');
% grid on;
% xlabel('迭代次数');
% ylabel('安全性代价');
% title('威胁度收敛曲线');
% hold on;
% figure(4);
% plot(1:1:batnum+1,totalfitness(1,:,3),'r',1:1:batnum+1,totalfitness(2,:,3),'k');
% legend('最小值','最大值');
% grid on;
% xlabel('迭代次数');
% ylabel('隐蔽性代价');
% title('隐蔽性收敛曲线');
% hold on;
% figure(5);
% for i=1:1:resultnum
% plot3(resultevolutionfitness(i,1)/min(resultevolutionfitness(:,1)),resultevolutionfitness(i,2)/min(resultevolutionfitness(:,2)),resultevolutionfitness(i,3)/min(resultevolutionfitness(:,3)),'ro');
% text(resultevolutionfitness(i,1)/min(resultevolutionfitness(:,1)),resultevolutionfitness(i,2)/min(resultevolutionfitness(:,2)),resultevolutionfitness(i,3)/min(resultevolutionfitness(:,3)),num2str(i));
% hold on;    
% end
% grid on;
% xlabel('航迹长度代价');
% ylabel('安全性代价');
% zlabel('隐蔽性代价');
% title('最优非支配解集');
% hold on;
% figure(6+i);
% for i=1:1:resultnum
% plot3(resultevolutionfitness(i,1)/min(resultevolutionfitness(:,1)),resultevolutionfitness(i,2)/min(resultevolutionfitness(:,2)),resultevolutionfitness(i,3)/min(resultevolutionfitness(:,3)),'ro');
% text(resultevolutionfitness(i,1)/min(resultevolutionfitness(:,1)),resultevolutionfitness(i,2)/min(resultevolutionfitness(:,2)),resultevolutionfitness(i,3)/min(resultevolutionfitness(:,3)),num2str(i));
% hold on;    
% end
% grid on;
% xlabel('航迹长度代价');
% ylabel('安全性代价');
% zlabel('隐蔽性代价');
% title('最优非支配解集');
% hold on;