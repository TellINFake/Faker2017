%TSP������Ⱥ�㷨

global NC;      %��������
global city_n;  %��������
global dis_table;   %���о������
global G;       %��¼��������
global everbest;    %�������Ž�
global hu_table;    %����ʽ�ֲ���
global adapt_best;

%=================================================================
%��һ��Ⱥ����
global a_A;       %��Ϣָ��
global b_A;       %����ʽָ��
global tobu_A;        %���Ͻ��ɱ�
global ph_table_A;    %��Ϣ�طֲ���
global ant_n_A;       %��������
global temp_pool_A;   %���ñ�ѡ����
global dispose_A;     %��Ϣ�ػӷ���
global best_individual_A; %�������Ž�·��
global adapt_ave_A;   %��������
global adapt_best_A;  %�������Ž�ı仯
global everbest_A;    %�������Ž�
%=================================================================
%�ڶ���Ⱥ����
global a_B;       %��Ϣָ��
global b_B;       %����ʽָ��
global tobu_B;        %���Ͻ��ɱ�
global ph_table_B;    %��Ϣ�طֲ���
global ant_n_B;       %��������
global temp_pool_B;   %���ñ�ѡ����
global dispose_B;     %��Ϣ�ػӷ���
global best_individual_B; %�������Ž�·��
global adapt_ave_B;   %��������
global adapt_best_B;  %�������Ž�ı仯
global everbest_B;    %�������Ž�
%=================================================================


initial;    %��ʼ�����ݣ�����������������λ��

for G=1:NC
    pause(0.01);
    search;         %�⹹��
    if mod(G,10)==0
        communication;  %ÿ10�ε�������һ��ͨ��
    end
    adapting;       %�������·����
    keepbest;       %���ź���
    ph_fresh;       %��Ϣ�ظ���
    paint;          %����ͼ��
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