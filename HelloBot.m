clear all
close all
clc

startup_rvc

cartesian = 0;

%% Definition of the links through the DH parameters

l1 = 0.147;
l2 = 0.155;
l3 = 0.135;
l4 = 0.2175;
offset = 0.033;
theta3_previous = 0;
L(1) = Link('d', l1,'a', offset, 'alpha', pi/2, 'offset', degtorad(169));
L(2) = Link('d',0,'a', l2, 'alpha', 0, 'offset', pi/2+degtorad(65));
L(3) = Link('d',0,'a', l3, 'alpha', 0, 'offset', -pi/2+degtorad(-146+90));
L(4) = Link('d', 0,'a', 0, 'alpha', pi/2, 'offset', pi);
L(5) = Link('d', l4,'a', 0, 'alpha', 0);
bot = SerialLink(L, 'name', 'HelloBot');

elbow_dir = -1;

%jointlimits
theta1Lim = [0.0100692 5.84014];  
theta2Lim = [0.0100692 2.61799];   
theta3Lim = [-5.02655 -0.015708];  
theta4Lim = [0.0221239 3.4297];    
theta5Lim = [0.110619 5.64159];

A = [0.0100692 0.0100692 -0.15708 0.0221239 0.110619];
T = bot.fkine(A);

%% Definition of the path
T = bot.fkine([0 0 0 0 0]);

% Points [x,y,z,theta,phi] to define the path, add in the c program the grip opening/closing
points = [0.008 -0.0016 0.4484 degtorad(99) -0.12;    %not insert it in the c program!!!!           
           0.15 0.1 0.3 degtorad(0) -3;
           0.3 0.2 0.4 degtorad(0) -3;
          -0.3 0.25 0.4 degtorad(0) -3;
          -0.3 0.3 0.2 degtorad(0) -3;
          -0.35 0.2 0.1 degtorad(0) -3;
          -0.4 0 0 degtorad(0) -3;
           -0.4 0 0 degtorad(0) -5.6;       
           -0.4 0 0 degtorad(0) -3;
           -0.35 0.2 0.1 degtorad(0) -3;
           -0.3 0.3 0.2 degtorad(0) -3;
           -0.3 0.25 0.4 degtorad(0) -3;
           0.3 0.2 0.4 degtorad(0) -3;
           ];
    
if (cartesian == 1)
    T = bot.fkine([0 0 0 0 0]);
    T2 = bot.fkine([degtorad(169), pi/2+degtorad(65), -pi/2+degtorad(-146+90), 0, 0]);
    
    %% Trajectory: to add point in the between of given points and have a better simulation
    N = 40;
    M = 0;
    P = 0;
    traj = zeros((size(points, 1) - 1)*N + 1, size(points, 2));
    q = zeros(size(traj, 1), size(traj, 2));
    q_degrees = zeros(size(traj, 1), size(traj, 2));
    
    for i = 1:(size(points, 1) - 1)    %i: number of given points
        for k = 0:(N-1)    %k: number of points to be adeed between two given points
            for j = 1:5    %j: number of degrees of freedom (without the gripper)
                traj(i + M, j) = points(i, j) + P*(points(i + 1, j) - points(i, j))/N;
            end
  
            if (k == (N-1))
            else
                M = M + 1;
            end
            P = P + 1;
            if (P > (N-1))
                P = 0;
            end
        end
    end
    
    traj(i + M + 1, :) = points(i + 1,:);
else
    traj = points;
end


%% Inverse dynamics
for i = 1:size(traj,1)
    % theta 1
    theta1 = atan2(traj(i, 2), traj(i, 1)) - degtorad(169) + pi;
    theta1 = -theta1;
    
    if (-theta1<=theta1Lim(1) || -theta1>=theta1Lim(2))
        disp(['theta1 out of bounds for i = ', num2str(i)])
    end
    
    t3 = sqrt(traj(i, 1)^2 + traj(i, 2)^2) - l4*cos(traj(i, 4)) - offset;
    z3 = traj(i, 3) - l4*sin(traj(i, 4)) - l1;
    
    % theta 3
    c3 = (t3^2 + z3^2 - l3^2 - l2^2)/(2*l3*l2);
    s3 = elbow_dir*sqrt(1-c3^2);
    theta3 = (atan2(s3, c3) + pi/2 - degtorad(-146+90));
    
    if (-theta3<=theta3Lim(1) || -theta3>=theta3Lim(2))
        disp(['theta3 out of bounds for i = ', num2str(i)])
    end
    
    % theta 2
    theta2 = atan2(z3, t3) - atan2(l3*s3, l2+l3*c3) - 0*pi/2 - degtorad(65) + 0*0.188;
    theta2 = (atan2(z3, t3) - atan2(l3*s3, l2+l3*c3) - pi/2 - degtorad(65) + 0.188);
 
    if (-theta2<=theta2Lim(1) || -theta2>=theta2Lim(2))
                disp(['theta2 out of bounds for i = ', num2str(i)])
    end
    
    % theta 4...
    theta4 = (traj(i, 4) - theta2 - theta3 - pi/2 - degtorad(9));

    if (-theta4<=theta4Lim(1) || -theta4>=theta4Lim(2))
        disp(['theta4 out of bounds for i = ', num2str(i)])
    end
    
    % theta 5..
    theta5 = traj(i, 5);
    if (-theta5<=theta5Lim(1) || -theta5>=theta5Lim(2))
            disp(['theta5 out of bounds for i = ', num2str(i)])
    end
    % put all the angles in a matrix
    q(i, :) = [theta1 theta2 theta3 theta4 theta5];

end
 
if (cartesian == 0)    % generation of a trajectory on the basis of the solution of the inverse kinematics. 
                       % (N) points are added to have a better simulation.
    N = 20;
    M = 0;
    P = 0;
    if (size(q,1) > 1)
        for i = 1:(size(q, 1) - 1)
            for k = 0:(N-1)    %k: number of points to be added betwenn two points
                for j = 1:5   %j: number of dof (no gripper) 
                    q_jtraj(i + M, j) = q(i, j) + P*(q(i + 1, j) - q(i, j))/N;
                end
                if (k == (N-1))
                else
                    M = M + 1;
                end
                P = P + 1;
                if (P > (N-1))
                    P = 0;
                end
            end
        end
        q_jtraj(i + M + 1, :) = q(i + 1,:);

                
        figure('units','normalized','outerposition',[0 0 1 1])
        bot.plot(q_jtraj, 'workspace', [-0.5 0.5 -0.5 0.5 -0.25 0.5]);
    else
        figure('units','normalized','outerposition',[0 0 1 1])
        bot.plot(q, 'workspace', [-0.5 0.5 -0.5 0.5 0 0.5]);
    end
end

q=-q;
q(1,:)=A;



% an out of boundaries for the first points is acceptable: they won't be
% implemented in the real robot, they are just used to start the
% simulation.
