% General Mesh Generator
% Author: Kanak Agarwal (ID25M802)
% Last Modified: 20th August, 2025
% Supports both equidistant and stretched meshes (both directions)
% Also supports different stretch factors in X and Y.

Lx = 1; % Length of the Domain in x
Ly = 1; % Length of the Domain in y

nx = 10; % Number of points in x
ny = 10; % Number of points in y

stretch_factor_X = 1.5;
stretch_factor_Y = 1.5;

% The summation of the cell widths i.e. dx's/dy's should be the length of 
% the domain (Lx/Ly). So given, the stretch factor, the number of cells 
% and the length of the domain it is possible to calculate the length of 
% the first cell from the summation formula of a geometric progression.
% These lengths are denoted by a_X and a_Y.

if stretch_factor_X == 1
    a_X = Lx/nx;
elseif stretch_factor_X <1
    a_X = Lx * (1 - stretch_factor_X)/(1 - stretch_factor_X^nx);
else
    a_X = Lx * (stretch_factor_X - 1)/(stretch_factor_X^nx - 1);
end

if stretch_factor_X == 1
    a_Y = Ly/ny;
elseif stretch_factor_X <1
    a_Y = Ly*(1 - stretch_factor_Y)/(1 - stretch_factor_Y^ny);
else
    a_Y = Ly*(stretch_factor_Y - 1)/(stretch_factor_Y^ny - 1);
end

dx = [a_X];
dy = [a_Y];

for i = 2:nx % from the second element/cell
    dx = [dx, (dx(i-1) * stretch_factor_X)]; % Array of dx's
end

for i = 2:ny % from the second element/cell
    dy = [dy, (dy(i-1) * stretch_factor_Y)]; % Array of dy's
end

grid_lines_X = [0];
grid_lines_Y = [0];

centre_X = [];
centre_Y = [];

x = 0;
y = 0;

for i = 1:length(dx)

    x = x + dx(i);
    y = y + dy(i);

    grid_lines_X = [grid_lines_X, x];
    grid_lines_Y = [grid_lines_Y, y];
    
    centre_X = [centre_X, (grid_lines_X(i) + grid_lines_X(i+1))/2];
    centre_Y = [centre_Y, (grid_lines_Y(i) + grid_lines_Y(i+1))/2];

end

% Plotting
[X,Y] = meshgrid(grid_lines_X,grid_lines_Y);
[x,y] = meshgrid(centre_X, centre_Y);
plot(X,grid_lines_Y)
hold on
plot(grid_lines_X,Y')
scatter(x,y, marker="o")
