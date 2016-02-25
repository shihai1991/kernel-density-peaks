% Labels = szy_cluster_dp(dist, k, isPlot) 或 Labels = szy_cluster_dp(dist, k)
% 基于密度峰值的聚类算法，dist是所有数据点的距离矩阵，k是要聚类的目标数，
% isPlot==true表示需要画图，isPlot==false表示不要画图，默认值是false.
% Labels是一个2xSize(dist, 1)大小的矩阵，每一列对应一个数据点，
% 第一行表示没有halo控制的聚类结果，第二行表示有halo控制的聚类结果。
% isKernal=true 用kernal函数
function Labels = cluster_dp(data, k, isPlot, isKernal)
ND = max(data(:,2));
NL = max(data(:,1));
if (NL>ND)
    ND = NL;   %% 确保 DN 取为第一二列最大值中的较大者，并将其作为数据点总数
end

%% data 第一个维度的长度，相当于文件的行数（即距离的总个数）   
 N = size(data,1);
 
%% 初始化为零 
 for i = 1:ND
   for j = 1:ND
     dist(i,j) = 0;
     c(i,j) = 0;
   end
 end
 
 %% 构建核函数

%% 利用 data 为 dist 数组赋值，注意输入只存了 0.5*DN(DN-1) 个值，这里将其补成了满矩阵  （这里不考虑对角线元素 ）
for i = 1:N-1
    ix = data(i,1);
    iy = data(i,2);
    for  j=i+1:N;
        jx = data(j,1);
        jy = data(j,2);
        if isKernal == true
            dist(i,j) = kernal(ix,iy,ix,iy)-2*kernal(ix,iy,jx,jy)+ kernal(jx,jy,jx,jy);
            dist(j,i) = dist(i,j);
        else
            c(i,j) = sqrt((ix-jx)^2+(iy-jy)^2);
            c(j,i) = sqrt((ix-jx)^2+(iy-jy)^2);
        end
    end
end

for i = 1:N
    if isKernal == true
        dist(i,i) = 0;
    else
        c(i,i) = 0;
    end 
end

if nargin == 2
    isPlot = false;
else
    isPlot = true;
end

if isKernal == true
    N = size(dist, 1) * (size(dist, 2) - 1) / 2;
    ND = size(dist, 1);
else
    N = size(c, 1) * (size(c, 2) - 1) / 2;
    ND = size(c, 1);
end


%%确定dc
percent = 2;
% fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);

position = round(N*percent/100);  %% round 是一个四舍五入函数  
if isKernal == true
    sda = sort(squareform(dist));
    dc = sda(position);
else
    sda = sort(squareform(c));
    dc = sda(position);
end
%% 计算局部密度 rho (利用 Gaussian 核)  
 fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);

%% 将每个数据点的 rho 值初始化为零  
rho = [];
for i = 1:ND
    rho(i) = 0.;
end

% Gaussian kernel
for i = 1:ND-1
    for j = i+1:ND
        if isKernal == true % 用kernal距离来计算rho
            rho(i) = rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
            rho(j) = rho(j)+exp(-(dist(j,i)/dc)*(dist(j,i)/dc));
        else % 用普通距离来计算rho
            rho(i) = rho(i)+exp(-(c(i,j)/dc)*(c(i,j)/dc));
            rho(j) = rho(j)+exp(-(c(j,i)/dc)*(c(j,i)/dc));
        end
    end
end

%% 先求矩阵列最大值，再求最大值，最后得到所有距离值中的最大值
if isKernal == true
    maxd=max(max(dist));
else
    maxd=max(max(c));
end

%% 将 rho 按降序排列，ordrho 保持序
[rho_sorted, ordrho]=sort(rho,'descend');

%% 处理 rho 值最大的数据点
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;

%% 生成 delta 和 nneigh 数组
if isKernal == true
    for ii=2:ND
        delta(ordrho(ii))=maxd;
        for jj=1:ii-1
            if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
                delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
                nneigh(ordrho(ii))=ordrho(jj);
            end
        end
    end
else
    for ii=2:ND
        delta(ordrho(ii))=maxd;
        for jj=1:ii-1
            if(c(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
                delta(ordrho(ii))=c(ordrho(ii),ordrho(jj));
                nneigh(ordrho(ii))=ordrho(jj);
            end
        end
    end
end

%% 生成 rho 值最大数据点的 delta 值
delta(ordrho(1))=max(delta(:));

%% 决策图
 disp('Generated file:DECISION GRAPH');
 disp('column 1:Density');
 disp('column 2:Delta');

 fid = fopen('DECISION_GRAPH', 'w');
 for i=1:ND
     fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
 end
 
%% 选择一个围住类中心的矩形
 disp('Select a rectangle enclosing cluster centers')

%% 每台计算机，句柄的根对象只有一个，就是屏幕，它的句柄总是 0

%% >> scrsz = get(0,'ScreenSize')

%% scrsz =

%            1           1        1280         800

% 1280 和 800 就是你设置的计算机的分辨率，scrsz(4) 就是 800，scrsz(3) 就是 1280
% scrsz = get(0,'ScreenSize');

%% 人为指定一个位置，感觉就没有那么 auto 了 :-)
% figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);

%% ind 和 gamma 在后面并没有用到
for i=1:ND
    ind(i)=i;
    gamma(i)=rho(i)*delta(i);
end

%% 利用 rho 和 delta 画出一个所谓的“决策图”
if isPlot
%     subplot(2,1,1)
    figure(1) % 决策树
    tt = plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
    title ('Decision Graph','FontSize',15.0)
    xlabel ('\rho')
    ylabel ('\delta')
end


% subplot(2,1,1)
% rect = getrect(1);
%% getrect 从图中用鼠标截取一个矩形区域， rect 中存放的是
%% 矩形左下角的坐标 (x,y) 以及所截矩形的宽度和高度
% rhomin=rect(1);
% deltamin=rect(2);%% 作者承认这是个 error，已由 4 改为 2 了!

%% 初始化 cluster 个数
% NCLUST=0;

NCLUST = k;

%% cl 为归属标志数组，cl(i)=j 表示第 i 号数据点归属于第 j 个 cluster
%% 先统一将 cl 初始化为 -1
for i=1:ND
    cl(i)=-1;
end

[B, Index] = sort(gamma, 'descend');

% cl是每个数据点的所属类别号
% icl是所有聚类中心的序号
icl = Index(1:k);
cl(Index(1:k)) = 1:k;

%% 在矩形区域内统计数据点（即聚类中心）的个数
% for i=1:ND
%   if ( (rho(i)>rhomin) && (delta(i)>deltamin))
%       % 修改上面的判断条件，根据降序排列后gamma的前几个点变化情况来决定类别。
%      NCLUST=NCLUST+1;
%      cl(i)=NCLUST;  %% 第 i 号数据点属于第 NCLUST 个 cluster

%       icl(NCLUST)=i; %% 逆映射,第 NCLUST 个 cluster 的中心为第 i 号数据点
%   end
% end

% fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);

% disp('Performing assignation')

%assignation  将其他数据点归类
for i=1:ND 
    if (cl(ordrho(i))==-1)
        cl(ordrho(i))=cl(nneigh(ordrho(i)));
    end
end
%% 由于是按照 rho 值从大到小的顺序遍历,循环结束后, cl 应该都变成正的值了. 

 

%% 处理光晕点，halo这段代码应该移到 if (NCLUST>1) 内去比较好吧
%halo
for i=1:ND
    halo(i)=cl(i);
end

% 获取每一个 cluster 中平均密度的一个界 bord_rho
if (NCLUST>1)

    % 初始化数组 bord_rho 为 0,每个 cluster 定义一个 bord_rho 值
    for i=1:NCLUST
        bord_rho(i)=0.;
    end
    
     % 获取每一个 cluster 中平均密度的一个界 bord_rho
    for i=1:ND-1
        for j=i+1:ND
             %% 距离足够小但不属于同一个 cluster 的 i 和 j
            if ((cl(i)~=cl(j))&&((isKernal==true&&(dist(i,j)<=dc))||((isKernal==false&&c(i,j)<=dc))))
                rho_aver=(rho(i)+rho(j))/2.;%% 取 i,j 两点的平均局部密度
                if (rho_aver>bord_rho(cl(i)))
                    bord_rho(cl(i))=rho_aver;
                end
                if (rho_aver>bord_rho(cl(j)))
                    bord_rho(cl(j))=rho_aver;
                end
            end
        end
    end

%% halo 值为 0 表示为 outlier
    for i=1:ND
        if (rho(i)<bord_rho(cl(i)))
            halo(i)=0;
        end
    end

%% 逐一处理每个 cluster
for i=1:NCLUST
    nc=0; %% 用于累计当前 cluster 中数据点的个数
    nh=0; %% 用于累计当前 cluster 中核心数据点的个数
    for j=1:ND
        if (cl(j)==i)
            nc=nc+1;
        end
        if (halo(j)==i)
            nh=nh+1;
        end
    end
         fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
end

if isPlot
    cmap=colormap;
    for i=1:NCLUST
        ic=int8((i*64.)/(NCLUST*1.));
%         subplot(3,1,1);
%         hold on
        figure(2) % 中心点
        plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
        hold on
    end
%     subplot(3,1,2);
    figure(3) 
    disp('Performing 2D nonclassical multidimensional scaling');
    plot(data(:,1),data(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on
    xlabel ('X');
    ylabel ('Y');
    for i=1:ND
        A(i,1)=0.;
        A(i,2)=0.;
    end
    
    for i=1:NCLUST
        nn=0;
        ic=int8((i*64.)/(NCLUST*1.));
        for j=1:ND
            if (cl(j)==i) %cl为簇的数据点  halo为核心数据点
                nn=nn+1;
                A(nn,1)=data(j,1);
                A(nn,2)=data(j,2);
            end
        end
%         hold on
        figure(4)
        grid on
        plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',4,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
% plot3(A(1:nn,1),sqrt(2).*A(1:nn,1).*A(1:nn,2),A(1:nn,2),'o','MarkerSize',4,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));        
hold on
        if isKernal == true
%             title ('density peak based on kernel function','FontSize',15.0);
        else
%             title ('density peak','FontSize',15.0);
        end
    end
end

%% 计算原始聚类图 add by hshi
% % subplot(3,1,3);
% plot(data(:,1),data(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
% title ('the original clusters','FontSize',15.0);
% xlabel ('X');
% ylabel ('Y');
% 
% for i=1:ND
%     A(i,1)=0.;
%     A(i,2)=0.;
% end
% 
for i=1:NCLUST
    nn=0;
    ic=int8((i*64.)/(NCLUST*1.));
    for j=1:ND
        if (data(j,3)==i)
            nn=nn+1;
            A(nn,1)=data(j,1);
            A(nn,2)=data(j,2);
        end
    end
%     hold on
    figure(5)
    plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',4,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
    hold on
end

%% =============================
for i=1:ND
    Labels(:, i) = [cl(i) halo(i)]';
end
end
