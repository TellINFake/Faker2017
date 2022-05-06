%TSP问题蚁群算法

global NC;      %迭代次数
global city_n;  %城市数量
global dis_table;   %城市距离矩阵
global G;       %记录进化代数
global everbest;    %历代最优解
global hu_table;    %启发式分布表
global adapt_best;

%=================================================================
%第一种群数据
global a_A;       %信息指数
global b_A;       %启发式指数
global tobu_A;        %蚂蚁禁忌表
global ph_table_A;    %信息素分布表
global ant_n_A;       %蚂蚁数量
global temp_pool_A;   %放置备选城市
global dispose_A;     %信息素挥发率
global best_individual_A; %历代最优解路径
global adapt_ave_A;   %进化数据
global adapt_best_A;  %历代最优解的变化
global everbest_A;    %历代最优解
%=================================================================
%第二种群数据
global a_B;       %信息指数
global b_B;       %启发式指数
global tobu_B;        %蚂蚁禁忌表
global ph_table_B;    %信息素分布表
global ant_n_B;       %蚂蚁数量
global temp_pool_B;   %放置备选城市
global dispose_B;     %信息素挥发率
global best_individual_B; %历代最优解路径
global adapt_ave_B;   %进化数据
global adapt_best_B;  %历代最优解的变化
global everbest_B;    %历代最优解
%=================================================================


initial;    %初始化数据，将蚂蚁随机放入城市位置

for G=1:NC
    pause(0.01);
    search;         %解构成
    if mod(G,10)==0
        communication;  %每10次迭代进行一次通信
    end
    adapting;       %计算各回路长度
    keepbest;       %保优函数
    ph_fresh;       %信息素更新
    paint;          %绘制图形
end

result;

clear a;
clear b;
clear a_A;
clear a_B;
clear b_A;
clear b_B;
clear NC;
clear h;
clear i;
clear j;
clear k;
clear n;
clear min_pos;
clear min_dis;
clear r;
clear tobu_A;
clear tobu_B;
clear dis_sum;
clear dispose_A;
clear dispose_B;
clear dis;
clear ant_n_A;
clear ant_n_B;
clear x1;
clear x2;
clear y1;
clear y2;
clear temp_pool_A;
clear temp_pool_B;
clear city;
clear city_n;
clear adapt_A;
clear adapt_B;
clear ada_sum;
clear ada_temp;
clear G;
clear ans;
clear everbest_A;
clear everbest_B;
clear dis_table;
clear best_individual;
clear best_individual_A;
clear best_individual_B;