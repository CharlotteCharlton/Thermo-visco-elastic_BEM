T = readmatrix("output.txt");

N0 = 12; M0= 24; % number of elements in theta (N) and phi (M) directions

Q = N0*M0;
vec = zeros(N0, 7);
p = 1;
theta = pi/(2*N0);
phi = pi/M0+2*pi/M0;

c = 1;

for j = 1:4:4*Q
    % cartesian displacement (or traction for rigid) components for left sphere
    u1 = T(j);
    u2 = T(j+1);
    u3 = T(j+2);
    u4 = T(j+3);
    
    % cartesian displacement (or traction for rigid) components for right sphere
    u12 = T(j + 3*Q);
    u22 = T(j+1+ 3*Q);
    u32 = T(j+2+ 3*Q);
    u42 = T(j+3+ 3*Q);
    
    % radial displacement u_r (or radial traction for rigid)
    vec(p, 1) = u1*sin(theta)*cos(phi) + u2*sin(theta)*sin(phi) + u3*cos(theta);
    vec(p, 4) = u12*sin(theta)*cos(phi) + u22*sin(theta)*sin(phi) + u32*cos(theta);
    
    % angular displacement u_theta (or traction t_theta for rigid)
    vec(p, 5) = u12*cos(theta)*cos(phi) + u22*cos(theta)*sin(phi) - u32*sin(theta);
    vec(p, 3) = u1*cos(theta)*cos(phi) + u2*cos(theta)*sin(phi) - u3*sin(theta);

    % temperature dif (or heat flux for rigid)
    vec(p, 7) = u4;

    vec(p, 2) = theta;
    vec(p, 6) = phi;
        
    phi = phi + 2*pi/M0;
        
    if c == M0
        theta = theta + pi/N0;
        phi = pi/M0+2*pi/M0;
        c = 0;
    end

    c = c + 1;       
    p = p + 1;
        
end

% indices corresponding to the azimuthal order of the plotted elements
in = 6;
in2 = 18;
in3 = 24;

% plot of radial displacement/traction at the surface of the left sphere
figure(1)
hold on
scatter(vec(in:M0:end, 2), abs(((vec(in:M0:end, 1)))), 75, 'filled', 'MarkerFaceColor',[0, 0.4470, 0.7410], 'MarkerEdgeColor',[0, 0.4470, 0.7410])
scatter(vec(in2:M0:end, 2), abs((vec(in2:M0:end, 1))), 75, "square", 'filled', 'MarkerFaceColor',[0.8500, 0.3250, 0.0980], 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980])
scatter(vec(in3:M0:end, 2), abs((vec(in3:M0:end, 1))), 75, "^",'filled', 'MarkerFaceColor',[0.9290, 0.6940, 0.1250], 'MarkerEdgeColor', [0.9290, 0.6940, 0.1250])

% plot of the temperature field/heat flux on the surface of the left sphere
figure(2)
hold on
scatter(vec(in:M0:end, 2), abs(((vec(in:M0:end, 7)))), 75, 'filled', 'MarkerFaceColor',[0, 0.4470, 0.7410], 'MarkerEdgeColor',[0, 0.4470, 0.7410])
scatter(vec(in2:M0:end, 2), abs((vec(in2:M0:end, 7))), 75, "square", 'filled', 'MarkerFaceColor',[0.8500, 0.3250, 0.0980], 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980])
scatter(vec(in3:M0:end, 2), abs((vec(in3:M0:end, 7))), 75, "^",'filled', 'MarkerFaceColor',[0.9290, 0.6940, 0.1250], 'MarkerEdgeColor', [0.9290, 0.6940, 0.1250])