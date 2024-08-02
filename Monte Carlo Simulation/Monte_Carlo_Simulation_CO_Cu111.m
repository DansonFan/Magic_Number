% -- By Dingxin Fan, July 2023

clc, clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the distribution of CO on the M(111) surface.
% Here we only consider Cu (M=Cu)
% We assume CO is only at the top position of Cu
% and that CO cannot simultaneously be on two adjacent Cu atoms (distance 2.55 Å)

% Input parameters
xgrid_0 = 24; % Number of Cu in x direction, must be a multiple of 6
ygrid_0 = 24; % Number of Cu in y direction, must be even
a = 3.6; % Lattice constant in Å, For Cu, a = 3.6
coverage = 0.26; % Desired coverage
divx = 2; % Number of divisions along the x-axis of the large grid
divy = 2; % Number of divisions along the y-axis of the large grid
timeout = 5; % Maximum allowed time to generate a single point in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create grid
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

% Output some input parameters
num_Cu = xgrid_0*ygrid_0; % Ignoring boundary conditions
num_Cu_per_unit = xgrid_0*ygrid_0/(divx*divy); % Ignoring boundary conditions
long_range_radius = sqrt(num_Cu*area_per_M/pi);
disp("# of Cu atom in total:");
disp(num_Cu);
disp("# of Cu atom per unit:");
disp(num_Cu_per_unit);
disp("System equivalent radius (Å):");
disp(long_range_radius);

% Create figure window
figure;
axis tight manual;
ax = gca;
% Set background color to sky blue
ax.Color = [21/255, 105/255, 224/255];  % RGB color values
hold on;

% Calculate the coordinates of division points on x and y axes
xDivisions = linspace(0, gridSize(2), divx+1);
yDivisions = linspace(0, gridSize(1), divy+1);

% Draw division lines on x-axis
for i = 2:divx
    line([xDivisions(i), xDivisions(i)], [0, gridSize(1)], 'Color', 'g', 'LineStyle', '--');
end

% Draw division lines on y-axis
for i = 2:divy
    line([0, gridSize(2)], [yDivisions(i), yDivisions(i)], 'Color', 'g', 'LineStyle', '--');
end

% Fix the position of the first point to (0.0, 0.0)
%x = 0.0;
%y = 0.0;
%grid(1, 1) = 1;
%points = [y, x];   

% Set up video output
videoFilename = 'point_animation.mp4';
fullPath = fullfile(pwd, videoFilename);
disp(['Video file save path: ', fullPath]);
video = VideoWriter(videoFilename, 'MPEG-4');
video.FrameRate = 10;  % Set frame rate
open(video);

% Start timing
tic;

% Set different random seeds
rng('shuffle');

good = 0;
% Randomly generate the first point
    while good == 0
        n = randi([0, nMax]);  % Randomly generate an integer in the range [0, nMax]
        m = randi([0, mMax]);  % Randomly generate an integer in the range [0, mMax]
        if mod(n, 2) == 0 && mod(m, 2) == 0 || mod(n, 2) == 1 && mod(m, 2) == 1
            % When n is even and m is even, or n is odd and m is odd
            x = n * xincrease_2;
            y = m * yincrease;
            grid(1,1) = 1;
            points = [y,x];
            good = 1;
        else
            % Otherwise, continue to the next iteration
            continue;
        end
    end

pointIndex = 2;  % Index variable
timeoutFlag = false;  % Timeout flag
lastUpdatedTime = toc;  % Time of last successful point generation and update
while pointIndex <= numPoints && ~timeoutFlag
    validPoint = false;
    startTime = toc;  % Get current time
    while ~validPoint
        % Randomly generate the position of the next point      
        n = randi([0, nMax]);  % Randomly generate an integer in the range [0, nMax]
        m = randi([0, mMax]);  % Randomly generate an integer in the range [0, mMax]
        
        if mod(n, 2) == 0 && mod(m, 2) == 0 || mod(n, 2) == 1 && mod(m, 2) == 1
            % When n is even and m is even, or n is odd and m is odd
            x = n * xincrease_2;
            y = m * yincrease;
        else
            % Otherwise, continue to the next iteration
            continue;
        end
        
        % Check distance between new point and existing points
        distances = sqrt(sum((points - [y, x]).^2, 2)); % Check This 🐕
        
        % Check if new point's coordinates meet the criteria
        if all(distances > 1) && x <= gridSize(2) && y <= gridSize(1)
            validPoint = true;
            grid(floor(y*10)+1, floor(x*10)+1) = 1;
            points = [points; y, x]; %#ok<AGROW>
            pointIndex = pointIndex + 1;  % Update index variable
            
            % Update plot
            plot(points(:, 2), points(:, 1), 'wo', 'MarkerSize', 6.5,'MarkerFaceColor', 'white');
            title(['Point ', num2str(pointIndex-1), ' of ', num2str(numPoints)]);
            axis equal;
            axis([0 gridSize(2) 0 gridSize(1)]);
            drawnow;
            
            % Write the current frame to video
            frame = getframe(gcf);
            writeVideo(video, frame);

            % Pause for a short time to observe the drawing of each point
            pause(0.001);
            
            % Update the time of the last successful point generation and update
            lastUpdatedTime = toc;
        end
        
        % Check if the time to generate the point exceeds the allowed time
        currentTime = toc;
        if currentTime - lastUpdatedTime > timeout
            timeoutFlag = true;  % Timeout, set the timeout flag to true
            break;  % Exit the inner loop
        end
    end
    
    % Check if the time to generate the point exceeds the allowed time
    currentTime = toc;
    if currentTime - lastUpdatedTime > timeout
        timeoutFlag = true;  % Timeout, set the timeout flag to true
        disp('Unable to achieve the desired coverage, the program will end directly.');
        break;  % Exit the outer loop
    end
end

% End timing
elapsedTime = toc;

% Close video
close(video);

% Save coordinates of each point
save('point_coordinates.mat', 'points');

% Display runtime
disp(['Program runtime: ', num2str(elapsedTime), ' seconds']);
disp(['Video saved as "', videoFilename, '"']);

% Calculate weight for each point
weights = ones(size(points, 1), 1); % Initialize all point weights to 1

% Ignoring boundary conditions
%{
% Adjust weights according to constraints
for i = 1:size(points, 1)
    if points(i, 1) == 0 || points(i, 2) == 0
        % If either x or y of the point is 0, weight is 0.5
        weights(i) = 0.5;
    end
    
    if points(i, 1) == 0 && points(i, 2) == 0
        % If both x and y are 0, weight is 0.25
        weights(i) = 0.25;
    end
end
%}
% Calculate total number of points
numPoints_stat = sum(weights);
numCO = numPoints_stat;
% Output total number of points
disp("# of CO in total:");
disp(numCO);
numPoints_stat = size(points, 1);  % Total number of points (before considering boundary conditions)
distancesPoints = zeros(numPoints_stat, 1);  % Store distances between each point and nearest point

% Calculate distances between points
for i = 1:numPoints_stat
    currentPoint = points(i, :);
    otherPoints = points([1:i-1, i+1:end], :);  % Other points excluding the current point

    % Calculate distances between the current point and other points
    pointDistances = sqrt(sum((otherPoints - currentPoint).^2, 2));

    % Find points within the threshold distance
    standard = minDistance+0.01;
    closePoints = otherPoints(pointDistances < standard, :);

    % Draw lines connecting points
    for j = 1:size(closePoints, 1)
        line([currentPoint(2), closePoints(j, 2)], [currentPoint(1), closePoints(j, 1)], 'Color', 'w','LineWidth', 1.8);
    end

    % Find the minimum distance
    minDistance_store = min(pointDistances);

    % Store the minimum distance
    distancesPoints(i) = minDistance_store;
end

% Count CO in each division region
regionCounts = zeros(divy, divx); % Initialize CO count per region to zero matrix
for i = 1:size(points, 1)
    xIndex = ceil(points(i, 2) / xDivisions(2));
    yIndex = ceil(points(i, 1) / yDivisions(2));
    
    % Check if index exceeds boundary, if so, set index to boundary value
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

% Output CO count per region
%disp("CO counts per region:");
%disp(regionCounts);

%num_Cu_circle = floor((ygrid*xgrid)/area_per_M);
%num_Cu = xgrid_0*ygrid_0-((xgrid_0-1)*0.5+(ygrid_0*0.5-1)*0.5)-0.75;
% Reorder rows
reorderedregionCounts = regionCounts(end:-1:1, :);
CO_Counts_store = reorderedregionCounts(:);
%overallCoverage = numCO/(num_Cu);
overallCoverage = numPoints/(xgrid_0*ygrid_0); % Ignoring boundary conditions
%regionCoverage = reorderedregionCounts/(num_Cu/(divx*divy));
regionCoverage = reorderedregionCounts/(xgrid_0*ygrid_0/(divx*divy));
Coverage_store = regionCoverage(:);

% Output statistical results
disp("Overall coverage:");
fprintf('%.2f%% ', overallCoverage * 100);
fprintf('\n');
disp("CO count per region:");
disp(reorderedregionCounts);
disp("Coverage per region:");
for i = 1:size(regionCoverage, 1)
    for j = 1:size(regionCoverage, 2)
        fprintf('%.2f%% ', regionCoverage(i, j) * 100);
    end
    fprintf('\n');
end

% Get current time
current_time = datestr(now);

% Output current time
disp(['Current time: ' current_time]);
