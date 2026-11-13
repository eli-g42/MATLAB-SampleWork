%% Description:
%{
Evaluates the unmodified Newtonian aerodynamics model, for
a given mesh in .dat format. Computes and plots the coefficient of pressure for each
surface element, and the coefficient of lift and drag for the entire geometry. The angle of
attack, freestream conditions, and grid are inputs to the code. 
%}

%% Reads the data and sorts into a 3D array for surface elements
clc; clear; close all;

%% Freestream properties
u_inf = 5100; % m/s
rho_inf = 3.0968e-4; % kg/m^3
p_inf = 21.955; % Pa

%% Geometry
fprintf('Select which geometry to analyze\n');
fprintf('(1) Flat Plate\n');
fprintf('(2) Sphere\n');
fprintf('(3) Apollo\n');

geometry = input('Enter your desired geometry: 1, 2, or 3\n');
 
switch geometry
    case 1
        grid_fname = 'Flat_Plate_Processed.dat';
        A_ref = 1; % m^2

    case 2
        grid_fname = 'Sphere_Processed.dat';
        r_sphere = 1; % m
        A_ref = pi*r_sphere^2; % m^2

    case 3
        grid_fname = 'Apollo_Processed.dat';
        r_apollo = 1.9558; % m
        A_ref = pi*r_apollo^2; % m^2

    otherwise
        error('Error! Please choose 1, 2, or 3\n')
end

%% Angle of Attack:
alpha_array = input('Enter an array of angles of attack in degrees (eg [0 20 40]\n');

if isempty(alpha_array) || ~isnumeric(alpha_array)
    error('Error! Please enter an array of angles of attack\n');
end

n_alpha = length(alpha_array);

%% ===== Read and Process Grid
% -- One row per node, first 3 nodes form 1st triangle, and so on
grid_raw = dlmread(grid_fname);

n_tri = size(grid_raw,1)/3;
if (mod(n_tri,1)~=0)
   error('Error, number of nodes is not divisible by 3\n');
end

% -- grid_x(k,j,i) has k-th dimension of j-th node of surf tri i
grid_x = zeros(3,3,n_tri);

for i=1:n_tri
   for j=1:3
      grid_x(1:3,j,i) = grid_raw( 3*(i-1)+j, 1:3);
   end
end

%% Calculate normal and triangle areas

normal_array = zeros(3, n_tri);
area_array = zeros(1, n_tri);

for i = 1:n_tri
    x1 = grid_x(:, 1, i);
    x2 = grid_x(:, 2, i);
    x3 = grid_x(:, 3, i);

    % calculating normal vector
    top = cross(x2-x1, x3-x1);
    bottom = norm(top);

    % inward normal
    if geometry == 1
        normal_array(:,i) = top/bottom;
    else
        normal_array(:, i) = -top/bottom;
    end

    area_array(i) = 0.5 * bottom;
end

%% Compute Exposure and Coefficient of Pressure for Each Surface Element

% Column = each angle of attack; Row = each element
exposed_array = false(n_tri, n_alpha);
cp_array = zeros(n_tri, n_alpha);
pressure_diff = zeros(n_tri, n_alpha);
A_per_element = zeros(n_tri, n_alpha);
N_per_element = zeros(n_tri, n_alpha);

A = zeros(1, n_alpha);
N = zeros(1, n_alpha);
L = zeros(1, n_alpha);
D = zeros(1, n_alpha);

c_L = zeros(1, n_alpha);
c_D = zeros(1, n_alpha);

for i = 1:n_alpha
    alpha = alpha_array(i);
    
    u_vec = u_inf * [cosd(alpha); -sind(alpha); 0];

    for j = 1:n_tri
        u_dot_normal = dot(u_vec, normal_array(:,j));

        if u_dot_normal > 0
            exposed_array(j,i) = true;

            cp_array(j,i) = 2*(u_dot_normal/u_inf)^2; % correct

            pressure_diff(j,i) = (0.5 * rho_inf * u_inf^2) * cp_array(j,i);

            A_per_element(j,i) = pressure_diff(j,i) * normal_array(1,j) * area_array(j);
            N_per_element(j,i) = pressure_diff(j,i) * normal_array(2,j) * area_array(j);
            

        end
    end

    A(i) = sum(A_per_element(:,i));
    N(i) = sum(N_per_element(:,i));

    L(i) = -N(i)*cosd(alpha) - A(i)*sind(alpha);
    D(i) = A(i)*cosd(alpha) - N(i)*sind(alpha);

    c_L(i) = L(i) / (0.5 * rho_inf * u_inf^2 * A_ref);
    c_D(i) = D(i) / (0.5 * rho_inf * u_inf^2 * A_ref);

    % Print results for current angle
    fprintf('\nAngle of attack: %.2f degrees\n',alpha);
    fprintf('Exposed triangles: %d\n',sum(exposed_array(:,i)));
    fprintf('Unexposed triangles: %d\n',sum(~exposed_array(:,i)));
    fprintf('Total triangles: %d\n', n_tri);
    fprintf('Nonzero Cp values: %d\n',sum(cp_array(:,i) ~= 0));
    fprintf('Zero Cp values: %d\n',sum(cp_array(:,i) == 0));

    fprintf('Coefficient of Lift: %.4f\n',c_L(i));
    fprintf('Coefficient of Drag: %.4f\n',c_D(i));

    % Sanity Check 2: Coefficient of Pressure Limits
    if any(cp_array(:) < 0) || any(cp_array(:) > 2)
        error('Error: Cp values are outside the theoretical limits.');
    end
    
    if any(cp_array(~exposed_array) ~= 0)
        error('Error: Unexposed elements have nonzero Cp values.');
    end
    
    fprintf('\nSanity Check 2 Passed: Cp values are within theoretical limits.\n\n');

end

% Sanity Check 1: Flat Plate at Perpendicular Flow
if geometry == 1 && any(alpha_array == 90)

    if any(abs(cp_array(:,alpha_array == 90)) ~= 2)
        error('Error: Flat plate Cp should equal 2 for perpendicular flow.');
    end

    fprintf('\nSanity Check 1 Passed: Flat plate Cp = 2 for perpendicular flow.\n');
end


%% Plotting c_p

% Triangle connectivity
tri_conn = zeros(n_tri,3);

for i = 1:n_tri
    tri_conn(i,:) = [3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3];
end

% Generate a plot for each angle of attack
for i = 1:n_alpha

    figure;

    trisurf(tri_conn, grid_raw(:,1), grid_raw(:,2), ...
        grid_raw(:,3), cp_array(:,i));

    title(sprintf('Coefficient of Pressure: Geometry %d (\\alpha = %.1f^\\circ)', ...
        geometry, alpha_array(i)));

    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');

    colorbar;
    clim([0 2]);
    axis equal;
    shading flat;
    view(3);

end

%% Part f
if geometry == 3
    fprintf('L/D: %d\n', L/D)
end