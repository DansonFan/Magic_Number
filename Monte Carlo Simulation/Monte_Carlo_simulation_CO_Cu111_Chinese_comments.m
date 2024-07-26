% -- By Dingxin Fan, July 2023

clc, clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%模拟CO在M(111)表面的分布。
%这里我们目前只考虑Cu (M=Cu)
%我们假设CO都在Cu的top位
%并且认为相邻的两个Cu原子上(距离2.55 Å)，无法同时站有CO

%Input parameters
xgrid_0 = 24; %x方向 # of Cu, 需要是6的倍数
ygrid_0 = 24; %y方向 # of Cu, 必须是偶数
a = 3.6; %lattice constant in Å, For Cu, a = 3.6
coverage = 0.26; % 想要达到的coverage
divx = 2; %要把大图x轴切分成的份数
divy = 2; %要把大图y轴切分成的份数
timeout = 5; %生成单个点允许的最大时间 in s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 创建网格
area_per_M = a^2*sqrt(3)/4;
minDistance = a*sqrt(6)/2;
xincrease = a*sqrt(2)/2;
xincrease_2 = 1.5*xincrease;
yincrease = a*sqrt(6)/4;
xlength = (xgrid_0-0.51)*xincrease;
xgrid = ceil(xlength);
ylength = (ygrid_0-0.9)*yincrease;
ygrid = ceil(ylength);
%xshift = xgrid/divx - xincrease;
%yshift = ygrid/divy - yincrease;
numPoints = ceil((xgrid_0*ygrid_0)*coverage);
gridSize = [ygrid, xgrid];
grid = zeros(gridSize);
nMax = (xgrid_0/6)*4-1;
mMax = ygrid_0;  

% 输出一些input的参数
num_Cu = xgrid_0*ygrid_0; %不考虑boundary的情况
num_Cu_per_unit = xgrid_0*ygrid_0/(divx*divy); %不考虑boundary的情况
long_range_radius = sqrt(num_Cu*area_per_M/pi);
disp("# of Cu atom in total:");
disp(num_Cu);
disp("# of Cu atom per unit:");
disp(num_Cu_per_unit);
disp("系统等价半径 (Å):");
disp(long_range_radius);

% 创建图像窗口
figure;
axis tight manual;
ax = gca;
% 设置背景色为天蓝色
ax.Color = [21/255, 105/255, 224/255];  % RGB颜色值
hold on;

% 计算x轴和y轴的分割点坐标
xDivisions = linspace(0, gridSize(2), divx+1);
yDivisions = linspace(0, gridSize(1), divy+1);

% 绘制x轴的分割线
for i = 2:divx
    line([xDivisions(i), xDivisions(i)], [0, gridSize(1)], 'Color', 'g', 'LineStyle', '--');
end

% 绘制y轴的分割线
for i = 2:divy
    line([0, gridSize(2)], [yDivisions(i), yDivisions(i)], 'Color', 'g', 'LineStyle', '--');
end

% 固定第一个点的位置为(0.0, 0.0)
%x = 0.0;
%y = 0.0;
%grid(1, 1) = 1;
%points = [y, x];   

% 设置视频输出
videoFilename = 'point_animation.mp4';
fullPath = fullfile(pwd, videoFilename);
disp(['视频文件保存路径：', fullPath]);
video = VideoWriter(videoFilename, 'MPEG-4');
video.FrameRate = 10;  % 设置帧率
open(video);

% 开始计时
tic;

% 设置不同的随机种子
rng('shuffle');

good = 0;
% 随机生成第一个点
    while good == 0
        n = randi([0, nMax]);  % 随机生成 [0, nMax] 范围内的整数
        m = randi([0, mMax]);  % 随机生成 [0, mMax] 范围内的整数
        if mod(n, 2) == 0 && mod(m, 2) == 0 || mod(n, 2) == 1 && mod(m, 2) == 1
            % 当 n 是负偶数且 m 是偶数，或者 n 是奇数且 m 是奇数时
            x = n * xincrease_2;
            y = m * yincrease;
            grid(1,1) = 1;
            points = [y,x];
            good = 1;
        else
            % 其他情况均作废，继续下一次循环
            continue;
        end
    end

pointIndex = 2;  % 索引变量
timeoutFlag = false;  % 超时标志
lastUpdatedTime = toc;  % 上一个点成功生成并更新的时间
while pointIndex <= numPoints && ~timeoutFlag
    validPoint = false;
    startTime = toc;  % 获取当前时间
    while ~validPoint
        % 随机生成下一个点的位置      
        n = randi([0, nMax]);  % 随机生成 [0, nMax] 范围内的整数
        m = randi([0, mMax]);  % 随机生成 [0, mMax] 范围内的整数
        
        if mod(n, 2) == 0 && mod(m, 2) == 0 || mod(n, 2) == 1 && mod(m, 2) == 1
            % 当 n 是负偶数且 m 是偶数，或者 n 是奇数且 m 是奇数时
            x = n * xincrease_2;
            y = m * yincrease;
        else
            % 其他情况均作废，继续下一次循环
            continue;
        end
        
        % 检查新点与已有点之间的距离
        distances = sqrt(sum((points - [y, x]).^2, 2)); %Check This 🐕
        
        % 检查新点的坐标是否满足条件
        if all(distances > 1) && x <= gridSize(2) && y <= gridSize(1)
            validPoint = true;
            grid(floor(y*10)+1, floor(x*10)+1) = 1;
            points = [points; y, x]; %#ok<AGROW>
            pointIndex = pointIndex + 1;  % 更新索引变量
            
            % 更新绘图
            plot(points(:, 2), points(:, 1), 'wo', 'MarkerSize', 6.5,'MarkerFaceColor', 'white');
            title(['Point ', num2str(pointIndex-1), ' of ', num2str(numPoints)]);
            axis equal;
            axis([0 gridSize(2) 0 gridSize(1)]);
            drawnow;
            
            % 将当前图像帧写入视频
            frame = getframe(gcf);
            writeVideo(video, frame);

            % 暂停一段时间，以便观察每个点的绘制过程
            pause(0.001);
            
            % 更新上一个点成功生成并更新的时间
            lastUpdatedTime = toc;
        end
        
        % 检查生成点的时间是否超过规定时间
        currentTime = toc;
        if currentTime - lastUpdatedTime > timeout
            timeoutFlag = true;  % 超时，将超时标志设置为true
            break;  % 跳出内循环
        end
    end
    
    % 检查生成点的时间是否超过规定时间
    currentTime = toc;
    if currentTime - lastUpdatedTime > timeout
        timeoutFlag = true;  % 超时，将超时标志设置为true
        disp('无法达到设定的coverage，程序将直接结束。');
        break;  % 跳出外循环
    end
end

% 结束计时
elapsedTime = toc;

% 关闭视频
close(video);

% 保存每个点的坐标
save('point_coordinates.mat', 'points');

% 显示运行时间
disp(['程序运行时间：', num2str(elapsedTime), ' 秒']);
disp(['视频已保存为 "', videoFilename, '"']);

% 计算每个点的权重
weights = ones(size(points, 1), 1); % 初始化所有点的权重为1

% 不考虑boundary的情况
%{
% 根据限制条件调整权重
for i = 1:size(points, 1)
    if points(i, 1) == 0 || points(i, 2) == 0
        % 如果点的x或者y两者之一等于0，权重为0.5
        weights(i) = 0.5;
    end
    
    if points(i, 1) == 0 && points(i, 2) == 0
        % 如果点的x和y同时等于0，权重为0.25
        weights(i) = 0.25;
    end
end
%}
% 计算总点数
numPoints_stat = sum(weights);
numCO = numPoints_stat;
% 输出总点数
disp("# of CO in total:");
disp(numCO);
numPoints_stat = size(points, 1);  % 总点数 (考虑边界条件之前)
distancesPoints = zeros(numPoints_stat, 1);  % 存储每个点与最近点的距离

% 计算点之间的距离
for i = 1:numPoints_stat
    currentPoint = points(i, :);
    otherPoints = points([1:i-1, i+1:end], :);  % 不包含当前点的其他点

    % 计算当前点与其他点之间的距离
    pointDistances = sqrt(sum((otherPoints - currentPoint).^2, 2));

    % 找到距离小于阈值的点
    standard = minDistance+0.01;
    closePoints = otherPoints(pointDistances < standard, :);

    % 绘制连线
    for j = 1:size(closePoints, 1)
        line([currentPoint(2), closePoints(j, 2)], [currentPoint(1), closePoints(j, 1)], 'Color', 'w','LineWidth', 1.8);
    end

    % 找到最小距离
    minDistance_store = min(pointDistances);

    % 存储最小距离
    distancesPoints(i) = minDistance_store;
end

% 根据分割区域统计CO数量
regionCounts = zeros(divy, divx); % 初始化每个区域的CO数量为零矩阵
for i = 1:size(points, 1)
    xIndex = ceil(points(i, 2) / xDivisions(2));
    yIndex = ceil(points(i, 1) / yDivisions(2));
    
    % 检查索引是否超出边界，如果超出边界，则将索引设置为边界值
    if xIndex > divx
        xIndex = divx;
    elseif xIndex < 1
        xIndex = 1;
    end
    
    if yIndex > divy
        yIndex = divy;
    elseif yIndex < 1
        yIndex = 1;
    end
    
    regionCounts(yIndex, xIndex) = regionCounts(yIndex, xIndex) + weights(i);
end

% 输出每个区域的CO数量
%disp("CO counts per region:");
%disp(regionCounts);

%num_Cu_circle = floor((ygrid*xgrid)/area_per_M);
%num_Cu = xgrid_0*ygrid_0-((xgrid_0-1)*0.5+(ygrid_0*0.5-1)*0.5)-0.75;
% 重新排列行
reorderedregionCounts = regionCounts(end:-1:1, :);
CO_Counts_store = reorderedregionCounts(:);
%overallCoverage = numCO/(num_Cu);
overallCoverage = numPoints/(xgrid_0*ygrid_0); %不考虑boundary的情况
%regionCoverage = reorderedregionCounts/(num_Cu/(divx*divy));
regionCoverage = reorderedregionCounts/(xgrid_0*ygrid_0/(divx*divy));
Coverage_store = regionCoverage(:);

% 输出统计结果
disp("总体coverage:");
fprintf('%.2f%% ', overallCoverage * 100);
fprintf('\n');
disp("每个区域内CO的数量:");
disp(reorderedregionCounts);
disp("每个区域内的coverage:");
for i = 1:size(regionCoverage, 1)
    for j = 1:size(regionCoverage, 2)
        fprintf('%.2f%% ', regionCoverage(i, j) * 100);
    end
    fprintf('\n');
end

% 获取当前时间
current_time = datestr(now);

% 输出当前时间
disp(['当前时间：' current_time]);
