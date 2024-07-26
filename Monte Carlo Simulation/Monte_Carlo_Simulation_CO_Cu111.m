clc, clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the distribution of CO on the M(111) surface.
% Currently, we only consider Cu (M=Cu)
% We assume that CO is only on the top sites of Cu
% And we assume that two adjacent Cu atoms (distance 2.55 Å) cannot both have CO

% Input parameters
xgrid_0 = 24; % Number of Cu atoms in x-direction, needs to be a multiple of 6
ygrid_0 = 24; % Number of Cu atoms in y-direction, must be even
a = 3.6; % Lattice constant in Å, For Cu, a = 3.6
coverage = 0.26; % Desired coverage
divx = 2; % Number of divisions along the x-axis
divy = 2; % Number of divisions along the y-axis
timeout = 5; % Maximum allowed time for generating a single point in seconds
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

% Calculate the coordinates of division points for x-axis and y-axis
xDivisions = linspace(0, gridSize(2), divx+1);
yDivisions = linspace(0, gridSize(1), divy+1);

% Draw x-axis division lines
for i = 2:divx
    line([xDivisions(i), xDivisions(i)], [0, gridSize(1)], 'Color', 'g', 'LineStyle', '--');
end

% Draw y-axis division lines
for i = 2:divy
    line([0, gridSize(2)], [yDivisions(i), yDivisions(i)], 'Color', 'g', 'LineStyle', '--');
end

% Fix the position of the first point at (0.0, 0.0)
%x = 0.0;
%y = 0.0;
%grid(1, 1) = 1;
%points = [y, x];   

% Set up video output
videoFilename = 'point_animation.mp4';
fullPath = fullfile(pwd, videoFilename);
disp(['Video file saved at: ', fullPath]);
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
            % Other cases are invalid, continue to the next iteration
            continue;
        end
    end

pointIndex = 2;  % Index variable
timeoutFlag = false;  % Timeout flag
lastUpdatedTime = toc;  % Time of the last successful point generation and update
while pointIndex <= numPoints && ~timeoutFlag
    validPoint = false;
    startTime = toc;  % Get the current time
    while ~validPoint
        % Randomly generate the next point's position      
        n = randi([0, nMax]);  % Randomly generate an integer in the range [0, nMax]
        m = randi([0, mMax]);  % Randomly generate an integer in the range [0, mMax]
        
        if mod(n, 2) == 0 && mod(m, 2) == 0 || mod(n, 2) == 1 && mod(m, 2) == 1
            % When n is even and m is even, or n is odd and m is odd
            x = n * xincrease_2;
            y = m * yincrease;
        else
            % Other cases are invalid, continue to the next iteration
            continue;
        end
        
        % Check the distance between the new point and existing points
        distances = sqrt(sum((points - [y, x]).^2, 2)); % Check This 🐕
        
        % Check if the new point's coordinates meet the conditions
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
        
        % Check if the time for generating the point exceeds the allowed time
        currentTime = toc;
        if currentTime - lastUpdatedTime > timeout
            timeoutFlag = true;  % Timeout, set timeout flag to true
            break;  % Break out of the inner loop
        end
    end
    
    % Check if the time for generating the point exceeds the allowed time
    currentTime = toc;
    if currentTime - lastUpdatedTime > timeout
        timeoutFlag = true;  % Timeout, set timeout flag to true
        disp('Unable to achieve the desired coverage. The program will terminate.');
        break;  % Break out of the outer loop
    end
end

% End timing
elapsedTime = toc;

% Close video
close(video);

% Save the coordinates of each point
save('point_coordinates.mat', 'points');

% Display runtime
disp(['Program runtime: ', num2str(elapsedTime), ' seconds']);
disp(['Video saved as "', videoFilename, '"']);

% Calculate the weight of each point
weights = ones(size(points, 1), 1); % Initialize weights of all points to 1

% Ignoring boundary conditions
%{
% Adjust weights according to boundary conditions
for i = 1:size(points, 1)
    if points(i, 1) == 0 || points(i, 2) == 0
        % If either x or y of the point is 0, set weight to 0.5
        weights(i) = 0.5;
    end
    
    if points(i, 1) == 0 && points(i, 2) == 0
        % If both x and y of the point are 0, set weight to 0.25
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
distancesPoints = zeros(numPoints_stat, 1);  % Store distance of each point
