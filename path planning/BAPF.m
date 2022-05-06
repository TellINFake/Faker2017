function [path, world, bestfitness3, iteration_fitness3]=BAPF(F3, size, radius, threat_count)
 % initialize obstacle world value
 load('XYZmesh.mat');
 world = worldInitialization('XYZmesh.mat');
 plotWorld(world);
 potential = SCinitialization();
 path = findMinimumPath(world,potential);
 plot(path(:,1),path(:,2));

%% worldInitialization
%% inilize values of obstacle world and start,end point 
function world = worldInitialization()

 % create obstacle world through initializing parameters
 initial_value.Size = 100;
 initial_value.Scorner = [0 0];
 initial_value.Ncorner = [initial_value.Size initial_value.Size];
 initial_value.NumObstacles = 10;
 world = createWorld(initial_value.NumObstacles,...
                     initial_value.Scorner,initial_value.Ncorner);
 % define the size of vehicle
 world.radius_vehicle = 1;                 
 % create start node and end node
 world.startNode = [0 0];
 node = generateRandEndnode(world);
 world.endNode = node;
 
%% SCinitialization
%% initialize the parameters of APF
function potential = SCinitialization()
 potential.kAttra = 5; % attractive force interia
 potential.kRepul = 150000;  % repulsive force interia
 potential.thresDistance = 5; % obstacle influence distance
 potential.step = 1;
 potential.vectorX = [1 0];
 potential.vectorY = [0 1]; 
 
%% createWorld
%% 创建障碍环境
function world = createWorld(NumObstacles,Scorner,Ncorner)

  % check to make sure that the region is nonempty
  if (Ncorner(1) <= Scorner(1)) || (Ncorner(2) <= Scorner(2))
      disp('Not valid corner specifications!');
      world=[];
      
  % create world data structure
  else
    world.NumObstacles = NumObstacles;
    world.Ncorner = Ncorner;
    world.Scorner = Scorner;
                          
    for i=1:NumObstacles,
        % create obstacle radius
        world.radius(i) = 5;
        % randomly pick center of obstacles
        cx = Scorner(1) + world.radius(i)...
            + (Ncorner(1)-Scorner(1)-2*world.radius(i))*rand;
        cy = Scorner(2) + world.radius(i)...
            + (Ncorner(2)-Scorner(2)-2*world.radius(i))*rand;
        world.cx(i) = cx;
        world.cy(i) = cy;
    end
  end

%% generate random node
function node = generateRandEndnode(world)
 % flag = 0 randendpoint
 % 如果将来想拓展产生其他的随机点 flag = 1时候就可以了
 % 随机endpoint的位置是在右上角的范围内产生的

 pe(1) = 4*world.Ncorner(1)/5 + (world.Ncorner(1) - 4*world.Ncorner(1)/5)*rand;
 pe(2) = 4*world.Ncorner(2)/5 + (world.Ncorner(2) - 4*world.Ncorner(2)/5)*rand;
 collisionState = collision(pe,0,world,[]);
 while (collisionState.flagRisk==1)
     pe(1) = 4*world.Ncorner(1)/5 ...
         + (world.Ncorner(1) - 4*world.Ncorner(1)/5)*rand;
     pe(2) = 4*world.Ncorner(2)/5 ...
         + (world.Ncorner(2) - 4*world.Ncorner(2)/5)*rand;
     collisionState = collision(pe,0,world,[]);
 end
 node = pe;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check if the input point belong to the unexpected region
function collisionState = collision(pe,flag,world,potential)
switch flag
    case 0   % endNode 的 collision判别
        % check if the position of end point satisfies the condition
        % condition one: if end point is on the top right corner
        % condition two: if end point and obstacles conflict
        % if above two conditions are satisfied, collision_flag = 0
        % namely, the random end point is well
        collisionState.flagRisk = 0;
        flag_center_distance = sqrt((world.cx - pe(1)).^2 ...
                                     + (world.cy - pe(2)).^2);
        flag_allowMinDistance = world.radius_vehicle + world.radius;
        flag_disComparation = flag_center_distance - flag_allowMinDistance;
        if ((pe(1)<3*world.Ncorner(1)/4) && (pe(2)<3*world.Ncorner(2)/4)...
                ||any(flag_disComparation<0))
            collisionState.flagRisk = 1;
        end
    case 1  % path过程中currentNode和障碍物collision的判别
        collisionState.flagRisk = 0;
        flag_center_distance = sqrt((world.cx - pe(1)).^2 ...
            + (world.cy - pe(2)).^2);
        thresCnDistance = potential.thresDistance + world.radius ...
            + world.radius_vehicle;
        if any(flag_center_distance<thresCnDistance)
            collisionState.flagRisk = 1;
            % mark the obstacle number risking vehicle
            numRisk = find(flag_center_distance<thresCnDistance);
            collisionState.numRisk = numRisk;
            collisionState.thresCnDistance = thresCnDistance(numRisk);
        end
end
 
 function path = findMinimumPath(world,potential)

 step = potential.step;
 currentPosition = world.startNode;
 Goal = world.endNode;
 Path(1,1:2) = currentPosition;
 i = 2;
 while(norm(currentPosition - Goal)>step) 
     forceRes = computeResultentForce(currentPosition,world,potential);
     Path(i,1) = currentPosition(1) + step*cosd(forceRes.angleAdvan);
     Path(i,2) = currentPosition(2) + step*sind(forceRes.angleAdvan);
     currentPosition = [Path(i,1),Path(i,2)];
     i = i+1;
 end
 path = Path;
     
%% compute_distance
%% the range of angleAdvan:[0 180] and [0 -180]
%% 求矩阵每一行的模：sqrt(dot(a,a，2))
%% compute the ditance between current position and obstacles
function forceRes = computeResultentForce(currentPosition,world,potential)
 
 Goal = world.endNode;
 vectorX = potential.vectorX;
 
 % 首先计算Attractive force 因为无论有没有斥力作用 只要还没到达目标就存在引力
 % vectorAttra：引力矢量相关参数：大小和方向
 vectorAttra = Goal - currentPosition; 
 angleAttra = acosd(dot(vectorAttra,vectorX)./norm(vectorAttra)/norm(vectorX));
 if vectorAttra(2)<0
     angleAttra = -angleAttra;
 end 
 forceAttra = potential.kAttra*norm(vectorAttra);
 forceAttraX = forceAttra*cosd(angleAttra);
 forceAttraY = forceAttra*sind(angleAttra);
 
 % 判断是否有障碍威胁
 collisionState = collision(currentPosition,1,world,potential);
 
 % 只有引力的作用 没有risk obstalces
 if collisionState.flagRisk == 0 
     forceRes.angleAdvan = angleAttra;
 
 %存在risk obstacles 有引力作用
 elseif collisionState.flagRisk == 1 
     numRisk = collisionState.numRisk;
     vectorX = repmat(vectorX,size(numRisk,2),1);
     % vectorRepul: 斥力矢量相关参数：障碍物到当前点矢量的大小和与X轴的方向
     vectorRepul(:,1) = (currentPosition(1) - world.cx(numRisk))';
     vectorRepul(:,2) = (currentPosition(2) - world.cy(numRisk))';
     distanceCnObstacle = sqrt(dot(vectorRepul,vectorRepul,2));
     angleRepul = acosd(dot(vectorRepul,vectorX,2)./distanceCnObstacle...
                      ./sqrt(dot(vectorX,vectorX,2)));
     angleRepul(vectorRepul(:,2)<0) = -angleRepul(vectorRepul(:,2)<0);
     thresCnDistance = collisionState.thresCnDistance';
     forceRepul = potential.kRepul...
                  .*(1./distanceCnObstacle-1./thresCnDistance)...
                  .*(1./(distanceCnObstacle.^2));
     forceRepulX = forceRepul.*cosd(angleRepul);
     forceRepulY = forceRepul.*sind(angleRepul);
     forceResX = sum(forceRepulX) + forceAttraX;
     forceResY = sum(forceRepulY) + forceAttraY;
     angleRes = atand(forceResY/forceResX);
     if forceResX < 0;
         angleRes = angleRes - 180;
     end
     forceRes.angleAdvan = angleRes;
 end

 iteration_fitness3=min(VectorX);

 if nargin<1,  para=[20 1000 0.5 0.5];  end
Fitness=0;
Sol=rand(0,1);
n=para(1);      % Population size, typically 10 to 40
N_gen=para(2);  % Number of generations of BA
A=para(3);      % Loudness  
r=para(4);      % Pulse rate
Qmin=0;         % Frequency minimum
Qmax=2;         % Frequency maximum
% Iteration parameters
N_iter=0;       % Total number of function evaluations
% Dimension of the search variables
d=10;           % Number of dimensions 
% Lower limit/bounds/ a vector
Lb=-2*ones(1,d);
% Upper limit/bounds/ a vector
Ub=2*ones(1,d);
% Initializing arrays
Q=zeros(n,1);   % Frequency
v=zeros(n,d);   % Velocities
itermax = 8000;            % maximum number of iterations (generations)
F = 0.6;                  % DE-stepsize F from interval [0, 2]
F1=0.45;
CR = 0.9;                 % crossover probability constant from interval [0, 1]
NP=10;
D=2;
XVmin=10;
XVmax=20;
w_min=0.1;
w_max=0.9;
u=4;
% Initialize the population/solutions
for i=1:n,
  Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
  Fitness(i)=Fun(Sol(i,:));
end
% Find the initial best solution
[fmin,I]=min(Fitness);
best=Sol(I,:);

% Start the iterations -- Bat Algorithm (essential part)  %
for t=1:N_gen, 
% Loop over all bats/solutions
        for i=1:n,
          P_success=Fitness(i)./i;
          w=w_min+(1-w_min./w_max)*P_success;
          Q(i)=Qmin+(Qmin-Qmax)*rand;
          v(i,:)=w*v(i,:)+(Sol(i,:)-best)*Q(i);
          S(i,:)=Sol(i,:)+v(i,:);
          % Apply simple bounds/limits
          Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
          Sol(i,:)=u*Sol(i-1,:)*(1-Sol(i-2,:))+Fitness(i-1)-iteration_fitness3;
          % Pulse rate
          if rand>r
          % The factor 0.001 limits the step sizes of random walks 
              S(i,:)=best+0.001*randn(1,d);
          end
     % Evaluate new solutions
       forceRes.angleAdvan=S(i,:);      
       Fnew=Fun(S(i,:));
     % Update if the solution improves, or not too loud
           if (Fnew<=Fitness(i)) && (rand<A) ,
                Sol(i,:)=S(i,:);
                Fitness(i)=Fnew;
           end

          % Update the current best solution
          if Fnew<=fmin,
                best=S(i,:);
                fmin=Fnew;
          end
        end
        N_iter=N_iter+n;
end
%-----DE Initialize--------------------------------------------------------
pop = zeros(NP,D);        % initialize pop
for i = 1:NP
   pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
end
popold = zeros(size(pop));% toggle population
val = zeros(1,NP);        % create and reset the "cost array"
bestmem = zeros(1,D);     % best population member ever
bestmemit = zeros(1,D);   % best population member in iteration
for i = 1:NP              % Evaluate the best member after initialization
    input = pop(i,:);
    output(i) = feval(fname,input);
end
val = output;
[bestvalit,idx] = min(val);
bestmemit = pop(idx,:);   % best member of current iteration
bestmem = bestmemit;      % best member ever
bestval = bestvalit;      % best value ever


tr = zeros(1,itermax);
for j = 1:itermax
    if (bestval-VTR) < 1e-5
        break
    end   
    popold = pop;
    % generate the trail population
    for i = 1:NP

        CR=0.1;
        if val(i)<=mean(val)
           CR=0.1+(0.6-0.1)*(val(i)-max(val))/(min(val)-max(val)); 
        end     

        % pick up the donor and differential vectors
        rd = fix(rand(1) * NP + 1);
        while val(rd) > val(i)
            rd = fix(rand(1) * NP + 1);
        end
        rb = fix(rand(1) * NP + 1);
        while rb == i || rb == rd
             rb = fix(rand(1) * NP + 1);
        end
        rc = fix(rand(1) * NP + 1);
        while rc == i || rc == rd || rc==rb
             rc = fix(rand(1) * NP + 1);
        end
        re = fix(rand(1) * NP + 1);
        while re == i || re == rd || re==rb  || re==rc
             re = fix(rand(1) * NP + 1);
        end
        % bulid a trial vector and crossover
        jr = fix(rand(1) * D + 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 找基变量
        tb=rd;
        if val(rb,:)<=min(val(rd,:),val(rc,:))
           tb=rb;
        end
        if val(rc,:)<=min(val(rd,:),val(rb,:))
           tb=rc;
        end
        % 找另两个差分向量;
        tm=rd;
        if val(rb,:)>=min(val(rd,:),val(rc,:)) && val(rb,:)<=max(val(rd,:),val(rc,:))
           tm=rb;
        end
        if val(rc,:)>=min(val(rd,:),val(rb,:)) && val(rc,:)<=max(val(rd,:),val(rb,:))
            tm=rc;
        end
        % 找最差的向量
        tw=rd;
        if val(rb,:)>=max(val(rd,:),val(rc,:))
           tw=rb;
        end
        if val(rc,:)>=max(val(rd,:),val(rb,:))
           tw=rc;
        end
        % 自适应缩放因子
        
        F=0.1+0.8*(val(tm)-val(tb))/(val(tw)-val(tb));

        for k = 1:D
            if rand(1)<=CR || k==jr
               diff = F * (popold(tm,k) - popold(tw,k));
               pop(i,k) = popold(tb,k) + diff;                
            else
                pop(i,k) = popold(i,k);
            end
        end
        for k = 1:D
           if pop(i,k) > XVmax(k)
               pop(i,k) = XVmax(k);
           end
           if pop(i,k) < XVmin(k)
               pop(i,k) = XVmin(k);
           end
        end
    end
end
function plotWorld(world)

 % plot obstacles
 for i=1:world.NumObstacles
     N = 3 + (7-3)*rand;
     th = 0:2*pi/N:2*pi;
     X = world.radius(i)*cos(th)*rand + world.cx(i);
     Y = world.radius(i)*sin(th)*rand + world.cy(i);
     world.handle_obstacle(i) = fill(X,Y,'k');
     hold on;
 end
 % plot start node and end node 
 N = 10; 
 th = 0:2*pi/N:2*pi;
 X = world.radius_vehicle*cos(th) + world.startNode(1);
 Y = world.radius_vehicle*sin(th) + world.startNode(2);
 fill(X,Y,'r');hold on;
 X = world.radius_vehicle*cos(th) + world.endNode(1);
 Y = world.radius_vehicle*sin(th) + world.endNode(2);
 fill(X,Y,'r');hold on;
 axis([0 100 0 100]);
 axis square;
 

 