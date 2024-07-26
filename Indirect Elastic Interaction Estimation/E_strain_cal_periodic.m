% -- Developed by Dingxin Fan, April 2024%

clc, clear

% Open the file
%fileID = fopen('path to your .xyz file','r');

% Skip the first two lines
for i = 1:2
    fgetl(fileID);
end

% Read the atomic information
atom_data = textscan(fileID, '%s %f %f %f');

% Parameters for CO/Cu(111)
mu_L = 0.00152; % Lame's second parameter (or shear modulus)
sigma = 0.34; % Poisson ratio
dilatation_2D = 0.1; %2D dilatation
xi = 1/3; % scaling factor
del_r = 3.51; % adsorption height
Hartree_to_eV = 27.211407953; %Convert Hartree to eV

% Close the file
fclose(fileID);

% Extract the data
element_names = atom_data{1};
x_coordinates = atom_data{2} * 1.8897259886; % Unit conversion (from Å to a.u.)
y_coordinates = atom_data{3} * 1.8897259886; % Unit conversion (from Å to a.u.)

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
num_cells_x = 2;
num_cells_y = 2;

% Define the CO separation threshold for interaction visualization
threshold = 10; % unit: a.u. (5.3nm = 100 a.u.)

% Plot the current unit cell and all surrounding unit cells
for k1 = -num_cells_x:num_cells_x
    for k2 = -num_cells_y:num_cells_y
        % Calculate the offset for neighboring unit cells
        delta_x = k1 * a1(1) + k2 * a2(1);
        delta_y = k1 * a1(2) + k2 * a2(2);
        
        % Plot the current unit cell
        scatter(C_x_coordinates + delta_x, C_y_coordinates + delta_y, 'filled');
        hold on;
        
        % Plot the boundaries of the current unit cell and its surrounding unit cells
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
% Note: the factor of 1/2 has been taken care of within for loops
E_elastic = 4 * pi * mu_L * (dilatation_2D)^2 * xi * ((del_r)^3)^2 * (1/(1-sigma)) * Hartree_to_eV * sum(rij_array .^ -3);
E_elastic_per_CO = E_elastic / num_C_atoms;

%disp(['E_elastic = ', num2str(E_elastic)]);
disp(['number of CO molecules = ', num2str(num_C_atoms)]);
disp(['E_elastic per CO (eV) = ', num2str(E_elastic_per_CO)]);

xlabel('X Coordinate');
ylabel('Y Coordinate');
title('CO Distribution and Interactions with Surrounding Unit Cells');
