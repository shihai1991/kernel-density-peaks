%% Initialization
clear ; close all; clc

% Load Training Data
fprintf('Loading Data\n');
% 399*3测试集
load('data1.mat');
% 3100*3测试集
% load('data2.mat'); 无法计算精确度
% 5300*3测试集
% load('data3.mat'); 无法计算精确度
% load('data4.mat');
% load('data5.mat');
% 鸢尾花数据集
% load('iris.mat');
% data1 = data(:,3);
% data2 = data(:,4);
% label = data(:,5);
% data = [data1, data2, label];
% 红酒数据集
% load('wine.mat');
% data1 = data(:,1);
% data2 = data(:,12);
% label = data(:,14);
% data = [data1, data2, label];
% glass数据集
% load('glass.mat');
% data1 = data(:,3);
% data2 = data(:,7);
% label = data(:,10);
% data = [data1, data2, label];

%% compute the data based on destiny peaks
% 聚类簇数
cluster_num = 5;
% cluster_num = max(data(:,3));
% 是否显示图集
isPlot = true;
% 是否使用核函数
 isKernal = true;
 isKernal = false;
% computing...
fprintf('computing...');
labels = cluster_dp(data, cluster_num, isPlot, isKernal); 

%% 对上述部分进行改进
% % 最大聚类簇数
% max_cluster_num = max(data(:,3));
% % 是否显示图集
% isPlot = true;
% % 是否使用核函数
%  isKernal = true;
% %  isKernal = false;
% % computing...
% fprintf('computing...');
% accury = zeros(1, max_cluster_num);
% time_compute = zeros(1, max_cluster_num);
% 
% for cluster_num = 2:max_cluster_num*2
%     time_before = now;
%     labels = cluster_dp(data, cluster_num, isPlot, isKernal);
%     time_after = now;
%     time_compute(cluster_num) = time_after-time_before;
%     accury(cluster_num) = AccMeasure(data(:,3), labels(1,:)');
% end
% figure(6)
% x = 2:cluster_num*2;
% plot(x,accury(x),'b--o');

%% compute the data based on k-means
% opts = statset('Display','final');
% [idx,C] = kmeans(data,cluster_num,'Distance','cityblock',...
%     'Replicates',5,'Options',opts);
% 
% figure;
% cmap = colormap;
% for i = 1:cluster_num
%     ic=int8((i*64.)/(cluster_num*1.));
%     hold on
%     plot(data(idx==i,1),data(idx==i,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
% end

%% 打印准确度
% fprintf('computing accuracy:%f', AccMeasure(data(:,3), labels(1,:)'));