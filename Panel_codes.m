%Script to output the lift coefficient for a given 4-digit NACA airfoil and 
%the specified flow conditions, as well as to plot the airflow around the 
%airfoil at a specified angle of attack
%
%The script asks for four user inputs to specified the NACA airfoil to use,
%the airspeed of the flow, the angle of attack of the airfoil and the
%number of panels to use to discretize the airfoil. The script then uses
%the two provided functions 'panelgen' and 'cdoublet'. 
%
%The script outputs the lift coefficient and a plot of the flow around the 
%airflow. The runtime of the script can be improved by reducting the 
%resolution of the domain by changing the X and Y values. On the contrary, 
%the resolution could be improved in the same way, with the sacrafice of 
%longer runtime. 




%clear the workspace to make sure there will be no unexpected problems in
%the script
clear
clc

%ask user for the flight parameters and number of panels to use
Airf = input('Which NACA airfoil to use? ','s');

%check if the NACA series is valid and raise an error if not
if any(isstrprop(Airf,'digit') ~= 1) || length(Airf) ~= 4
    error('Invalid NACA airfoil')
elseif isempty(Airf)
    error('Invalid NACA airfoil')
end

U_inf = input('What is the freestream velocity (m/s)? ','s');

%check if the velovity input is a number
if any(isstrprop(U_inf,'alpha')) == 1 || isempty(U_inf)
    error('Invalid velocity input - please provide a velocity in m/s')
else
    U_inf = str2double(U_inf);
end

AoA = input('What is the angle of attack (degrees)? ','s');

%check if the provided angle of attack is a number
if any(isstrprop(AoA,'alpha')) == 1 || isempty(AoA)
    error('Invalid angle input - please provide an angle in degrees')
elseif abs(str2double(AoA)) > 10
    warning('Large angle of attack provided - the output is likely to be inaccurate')
else
    AoA = str2double(AoA);
end

n = input('Number of panels to use: ');

%check if the no. of panels input is an integer
if mod(n,1) ~= 0 || n < 1 || isempty(n)
    error('Invalid panels input - please provide a positive, integer number')
end

tic
%convert the provided angle of attack into radians
AoA = deg2rad(AoA);

%use panelgen function to discretize the airfoil
[x,z]=panelgen(Airf,n,AoA); 

%calculate angle beta according to the equations provided in the handout
beta(1:n+1) = atan((z(2:n+2)-z(1:n+1))./(x(2:n+2)-x(1:n+1)));

%change the angles beta to a desired domain if necessary. Since atan()
%provides outputs between -pi/2 and pi/2, if angle is less than 0 we need
%to add 2pi
if beta < 0
    beta = beta + 2*pi;
end

%initialize the array 'midpoints', which will contain all midpoints of the
%created panels
midpoints = zeros(n+1,2);

%fill the array 'midpoints' with the x-coordinates of midpoints in the
%first column and the z-coordinates in the second column
midpoints(:,1) = (x(1:n+1)+x(2:n+2))/2;
midpoints(:,2) = (z(1:n+1)+z(2:n+2))/2;

%create a grid of arrays x and z, with the arrays being repeated over n+1
%rows. This will make it possible to use cdoublet with array opperations
x_grid = meshgrid(x,1:n+1);
z_grid = meshgrid(z,1:n+1);

%initialize arrays 'u' and 'v'
u = zeros(n+1,n+1);
v = zeros(n+1,n+1);

%iterate over columns in x_grid and z_grid arrays
for j = 1:n+1
    %each iteration of the loop takes all points from 'midpoints' array and
    %calls a modified 'cdoublet' function to calculate the velocities
    %induced on all points by a panel with given x/z endpoint coordinates, 
    %using array opperations
    [u(:,j),v(:,j)] = cdoublet_array([midpoints(:,1),midpoints(:,2)],[x_grid(:,j),z_grid(:,j)],[x_grid(:,j+1),z_grid(:,j+1)]); 
end

%create a length(beta)xlength(beta) grid of array beta, with the array 
%being repeated over n+1 rows, so as to allow array opperations
beta_grid = meshgrid(beta);

%transpose 'the beta_grid' matrix to fit the desired conditions
beta_grid = beta_grid';

%apply the equation provided in the handout to create a matrix that will
%form the equations for solving a linear system
A = v.*cos(beta_grid) - u.*sin(beta_grid);

%initialize the vector B
B = zeros(n+1,1);

%fill the vector B with constants for solving the system of simultaneous
%equations 
B([1:n],1) = -U_inf*sin(AoA-beta(1:n));

%implement the Kutta condition by creating the last equation and equating
%it to zero - done when initializing vector B
A(n+1,[1,n,n+1]) = [1,-1,1];

%initialize the solutions array 'mu'
mu = zeros(n+1,1);

%solve the system for n+1 variables 'mu' using the 'pinv()' function,
%equating the result to an array 'mu'
mu = linsolve(A,B);

%calculate the lift coefficient according to the equation provided in the
%handout
Cl = -2*mu(n+1)/U_inf;
%print the result of the calculation on screen
fprintf('CL = %f\n',Cl)


%define the domain and range over which the streamlines will be calculated
%the model can be improved by reducting the increment of arrays X and Y,
%but 0.01 has been sufficient when testing this script
X= -0.2:0.01:1.2;
Y = -0.7:0.01:0.7;

%initialize array 'points' making it length(X)^2 long, corresponding to the
%number of possible coordinates in the created domain (X/Y)
points = zeros(length(X)^2,2);

%create a grid of arrays X and Y
[X,Y] = meshgrid(X,Y);

%make the first and second column of array 'points' equal to all possible X
%and Y coordinates, respectively - each possible X is paired with each
%possible Y
points(:,1) = X(1:end);
points(:,2) = Y(1:end);

%define varialbe range to be equal to the number of rows of array 'points'
range = size(points,1);

%create a grid of arrays x and z, with the arrays being repeated over the
%number of rows of the array 'points'
x_grid_2 = meshgrid(x,1:range); 
z_grid_2 = meshgrid(z,1:range);

%initialize arrays 'u_all' and 'v_all'
u_all = zeros(range,n+1);
v_all = zeros(range,n+1);

%iterate over columns in x_grid_2 and z_grid_2 arrays
for q = 1:n+1
    %each iteration of the loop takes all points from 'points' array and
    %calls a modified 'cdoublet' function to calculate the velocities
    %induced on all points in desired domain by a panel with given x/z 
    %endpoint coordinates, using array opperations
    [u_all(:,q),v_all(:,q)] = cdoublet_array([points(:,1),points(:,2)],[x_grid_2(:,q),z_grid_2(:,q)],[x_grid_2(:,q+1),z_grid_2(:,q+1)]); % p is the point and p1,p2 are two end points
end

%create a grid of array 'mu_grid', with the arrays being repeated over the 
%number of rows of the array 'points'
mu_grid = meshgrid(mu,1:range);

%initialize arrays U and V
U = zeros(size(mu_grid));
V = zeros(size(mu_grid));

%calculate entries in arrays U and V as a product of mu and the velocity
%induced by a given panel
U = mu_grid.*u_all;
V = mu_grid.*v_all;

%sum the rows of the arrays U and V
U = sum(U,2);
V = sum(V,2);

%add a constant term to each entry in arrays U and V, according to equation
%1 from the handout
U = U+U_inf*cos(AoA);
V = V+U_inf*sin(AoA);

%split the arrays U and V so that they form a matrix with the number of
%columns equal to length(X)
U = reshape(U,length(X),[]);
V = reshape(V,length(X),[]);

%use function 'inpolygon()' to check if a given point at coordinate 
%(points(:,1),points(:,2)) is inside or on the polygon defined by arrays x
%and z
in = inpolygon(points(:,1),points(:,2),x,z);

%use the output of the function above to change the velocities at those
%coordinates to be equal to 'NaN', making them invisible on the plot
U(in==1) = nan;
V(in==1) = nan;

figure;
%plot the final velocities as streamlines using the 'streamslice()'
%function
subplot(1,2,1)
streamslice(X,Y,U,V)
hold on                     %hold on to be able to plot on the same graph
plot(x([1:n,1]),z([1:n,1]),'r')     %plot the airfoil using arrays x and z
hold off                    %finish plotting

subplot(1,2,2)
streamslice(X,Y,U,V,2)
hold on                     %hold on to be able to plot on the same graph
plot(x(1:n),z(1:n),'r')     %plot the airfoil using arrays x and z
hold off                   
%end of script