function [best,fmin,N_iter]=IBA(para)
% Display help
help bat_algorithm.m
%___________________
%Default parameters
if nargin<1, para=[20 100 0.5 0.5]; end
n=para(1);     % Population size,typically 10 to 40
iter_max=200;
N_gen=para(2); % Number of generations
A=para(3);     % Loundness (constant or decreasing)
r=para(4);     % Pulse rate (constant or decreasing)
% This frequency range determines the scalings
% You should change these values if necessary
Qmin=0;        % Frequency minimum
Qmax=2;        % Frequency maximum
% Iteration parameters
N_iter=0;      % Total number of function evaluations
% Dimension of the search variables
d=50;         % Number of dimensions
r0=r;
A0=A;
mean_max=0.95; % 随机权重平均的最大值
mean_min=0.5;  % 随机权重平均的最小值
sigma=0.1;     % 随机权重的方差
Low=-100;
High=100;
% Lower limit /bounds /a vector
Lb=Low*ones(1,d);
% Upper limit /bounds /a vector
Ub=High*ones(1,d);
% Initializing arrays
Q1=zeros(n,1); % Frequency one
Q2=zeros(n,1); % Frequency one
v=zeros(n,d);  % Velocities
% Initialize the population /solutions
Location=0;
for i=1:n
    Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
    Fitness(i)=Fun(Sol(i,:),d);
    Location=Location+Sol(i,:);
    % disp(['初始：',num2str(Sol(i,:)),'初始函数：'num2str(Fitness(i))]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始没问题%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%新添加限定%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
% Find the initial best solution
[fmin,I]=min(Fitness); %矩阵
avg_Location=Location/n;
avag=[];
avag=mean(Fitness);
best=Sol(I,:);
temp=[];

for t=1:N_gen
    %Loop over all bats/solutions
    sum=0; %计算适应值总和
    location=0;
    for i=1:n
        Q1(i)=1-exp(-abs(avag-fmin))*rand;
        Q2(i)=1-Q1(i);
        miu=mean_min+(mean_max-mean_min)./(iter_max)*rand(); %adaptive interia weight 
        w=miu+sigma*randn(); %location interia weight updating 
        tem(i,:)=Sol(i,:).*rand(1,d)*2; %generate a solution randomly
        %S(i,:)=w*Sol(i,:)+(avg_Location-Sol(i,:))*Q1(i)+(best-Sol(i,:))*Q2(i);%
        S(i,:)=w*Sol(i,:)+(tem(i,:)-Sol(i,:))*Q1(i)+(best-Sol(i,:))*Q2(i);
        % S(i,:)=Sol(i,:)+v(i,:);
        % Apply simple bounds/limits
        Sol(i,:)=simplebounds(S(i,:),Lb,Ub);
        location=location+Sol(i,:);
        % Pulse rate
        if rand>r0
            % The factor 0.001 limits the step sizes of random walks
            S(i,:)=best+A*randn(1,d);
        end
        S(i,:)=simplebounds(S(i,:),Lb,Ub);
        % Evaluate new solutions
        Fnew=Fun(S(i,:),d);
        % Update if the solution improves, or not too loud
        if (Fnew<=Fitness(i)&&(rand<A0))
            Sol(i,:)=S(i,:);
            Fitness(i)=Fnew;
            % 要增大ri, 减少Ai
        end
        % Update the current best solution
        if Fnew<=fmin
            best=S(i,:);%这个是解
            fmin=Fnew;
        end
        Sum=sum+Fnew;
        %disp(['位置：'，num2str(Sol(i,:)),'函数：'，num2str(Fnew)]);
    end
    avg_Location=Location/n;
    avag=sum/n;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%变化的%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A0=A*(exp(-0.1*t))*0.98/0.9;
    r0=r*(1-A0/A);
    % A0=A*0.9;
    % r0=r*(1-exp(-0.1*t));
    % disp(['A=',num2str(A0),'R='num2str(r0)]);
    N_iter=N_iter+n;
    temp=[temp, fmin];
    % plot(1:N_gen,fmin,'.');
end
plot(1:N_gen, log10(temp),'-*r');
xlabel('Times'); %标注横坐标
ylabel('log10(temp)'); %标注纵坐标
legend('FSABA');
hold on % 允许在同一坐标系下绘制不同的图形
% Output/display
disp(['Number of evaluations:', num2str(N_iter)]);
disp(['位置：',num2str(Sol(i,:)),'函数：',num2str(Fun(best,d))]);
% disp(['Best=',num2str(best),'fmin=',num2str(fmin)]);
% plot(1:N_gen,fmin,'-');

% Applicaiton of simple limits/bounds
function s=simplebounds(s,Lb,Ub)
% Apply the lower bound vector
for i=1:length(s)
    ns_tmp=s(:,i);
    if ns_tmp<Lb(:,i)
        ns_tmp=Lb(:,i)*rand;
    end
% Apply the upper bound vector
    if ns_tmp>Ub(:,i)
        ns_tmp=Ub(:,i);
    end
    %Update this new move
    s(:,i)=ns_tmp;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function: you own objective funciton can be written here
% Note: When you use your own function,please remember to change limits Lb
% and Ub (see lines52 to 55）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function z=Fun(u,d)

% 函数 F1 Sphere function with fmin=0 at(0,0,...,0)
function z=Fun(u,d)
z=sum(u.^2);
% 