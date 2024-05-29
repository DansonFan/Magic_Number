% -- Developed by Dingxin Fan, April 2024% -- Developed by Dingxin Fan, April 2024

clc, clear

% Open the file
fileID = fopen('path','r');
% Skip the first two lines
for i = 1:2
    fgetl(fileID);
end

% Read the atomic information
atom_data = textscan(fileID, '%s %f %f %f');

% Close the file
fclose(fileID);

% Extract the data
element_names = atom_data{1};
x_coordinates = atom_data{2} * 1.8897259886; % 单位转换
y_coordinates = atom_data{3} * 1.8897259886; % 单位转换

% Extract x and y coordinates of rows where the element is C
C_indices = strcmp(element_names, 'C');
C_x_coordinates = x_coordinates(C_indices);
C_y_coordinates = y_coordinates(C_indices);

% Get the number of C atoms
num_C_atoms = sum(C_indices);

% Add the distances between CO molecules in neighboring unit cells to the current unit cell
% Assume a 2D lattice, and "infinite" number of neighboring unit cells
a1 = [57.6773796603, 0];
a2 = [0, 49.9500760085];

% # of periodic unit cell in x and y directions
num_cells_x = 0;
num_cells_y = 0;

% Define the CO separation threshold for interaction visualization
threshold = 10; % unit: a.u. (5.3nm = 100 a.u.)

% 绘制当前 unit cell 及其周围的所有 unit cells
for k1 = -num_cells_x:num_cells_x
    for k2 = -num_cells_y:num_cells_y
        % 计算相邻 unit cell 的偏移量
        delta_x = k1 * a1(1) + k2 * a2(1);
        delta_y = k1 * a1(2) + k2 * a2(2);
        
        % 绘制当前 unit cell
        scatter(C_x_coordinates + delta_x, C_y_coordinates + delta_y, 'filled');
        hold on;
        
        % 绘制当前 unit cell 与其周围的所有 unit cells 的边界
        x_boundary = [min(C_x_coordinates) + delta_x, max(C_x_coordinates) + delta_x];
        y_boundary = [min(C_y_coordinates) + delta_y, max(C_y_coordinates) + delta_y];
        plot([x_boundary(1), x_boundary(2), x_boundary(2), x_boundary(1), x_boundary(1)], ...
             [y_boundary(1), y_boundary(1), y_boundary(2), y_boundary(2), y_boundary(1)], 'k--');
    end
end

% Calculate the distance between each pair of C atoms and name them
rij_array = [];
for i = 1:num_C_atoms
    for j = (i + 1):num_C_atoms
        dx = C_x_coordinates(j) - C_x_coordinates(i);
        dy = C_y_coordinates(j) - C_y_coordinates(i);
        rij = sqrt(dx^2 + dy^2);
        rij_array = [rij_array; rij];
        if rij < threshold
            plot([C_x_coordinates(i), C_x_coordinates(j)], [C_y_coordinates(i), C_y_coordinates(j)], 'r--');
        end
    end
end

for i = 1:num_C_atoms
    for j = 1:num_C_atoms
        for k1 = -num_cells_x:num_cells_x % Iterate over the offsets of neighboring unit cells
            for k2 = -num_cells_y:num_cells_y
                if k1 == 0 && k2 == 0
                    continue; % Skip the current unit cell
                end
                delta_x = k1 * a1(1) + k2 * a2(1);
                delta_y = k1 * a1(2) + k2 * a2(2);
                dx = C_x_coordinates(j) + delta_x - C_x_coordinates(i);
                dy = C_y_coordinates(j) + delta_y - C_y_coordinates(i);
                rij = sqrt(dx^2 + dy^2);
                rij_array = [rij_array; rij];
                if rij < threshold
                    plot([C_x_coordinates(i), C_x_coordinates(j)], [C_y_coordinates(i), C_y_coordinates(j)], 'r--');
                end
            end
        end
    end
end

% Calculate E_elastic per CO
E_elastic = 4.67 * pi * sum(rij_array .^ -3) / num_C_atoms;

disp(['E_elastic per CO = ', num2str(E_elastic)]);

xlabel('X Coordinate');
ylabel('Y Coordinate');
title('CO Distribution and Interactions with Surrounding Unit Cells');