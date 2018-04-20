clc; clear all; close all;

%Units across the work:
%Angle: radians until mentioned otherwise
%distance: milimeter
%Velocity: mm/sec
%Angular velocity: rad/sec
%Force: kg mm/sec^2

%% Question 1 Frames & DH Parameters
%Equating the parameters found 
%Symbolic joint angles
syms b c D e f g;
syms t1 t2 t3;
syms dt1 dt2 dt3;

DH = [b+D,        0, -c,     -pi/2;
    0, -pi/2+t1,  e,     0;
    0,       t2,  f,     0;
    0, -pi/2+t3,  0, -pi/2;
    g,        0,  0,     0];


N = size(DH, 1);
%% Question 2 Position Kinematics
%Generating Transformation matrices for each link
% Creating symbolic variables
dz =    DH(:,1);
theta = DH(:,2);
a =     DH(:,3);
alpha = DH(:,4);

% Calculating the intermediate transformations
T = sym(zeros(4,4,N));
for i = 1:N
    T(:,:,i) = dhparam2matrix(dz(i), theta(i), a(i), alpha(i));
    disp( T(:,:,i));
end

% This is the complete tranformation matrix from Fr to Ft
T04 = sym(eye(4,4));
for i = 1:N
    T04 = T04 * T(:,:,i);
end

disp("Final T matrix");
disp(T04);

%% Question 3 Velocity Kinematics
%Extracting xt from the final transformation matrix
xt = T04(1:3,4);

%Calculating Jacobian for translation 
% Jv = [diff(xt(1),t1), diff(xt(1),t2), diff(xt(1),t3); ...
%       diff(xt(2),t1), diff(xt(2),t2), diff(xt(2),t3); ...
%       diff(xt(3),t1), diff(xt(3),t2), diff(xt(3),t3)];
Jv = jacobian(xt, [t1,t2,t3]);

%For rotational jacobian
k = [0;0;1]; %standard
p = 1; %for revolute joint
r1_z0 = p*k; %z(i-1) = Rotmat of 0 to i-1
r2_z1 = p*(T(1:3,1:3,1)*T(1:3,1:3,2))*k; %R(0->1)
r3_z2 = p*(T(1:3,1:3,1)*T(1:3,1:3,2)*T(1:3,1:3,3))*k; % R(0->2)
Jw = [r1_z0, r2_z1, r3_z2];

J = [Jv;Jw];
disp("The 6 DOF Jacobian is:");
disp(J);

disp("3 DOF Jacobian:");
disp(Jv);

%As seen in the output below the second column is of the Jacobian is
%zero which means there's no movement in y-direction

%% Question 4 Force Propagation
%The torque tau = J' * F
disp("Symbolic relationship between the task space Wrench:")
sym_rel = J'

%% Question 5 Numeric Solution
b=424 ;c=300 ;D=380 ;e=328 ;f=323 ;g=82.4 ; %mm
t1=deg2rad(20);t2=deg2rad(90);t3=deg2rad(30); %rad
%5.a
disp("T matrix for the given configuration ");
out1 = double(subs(T04))

%5.b
joint_vel1 = degtorad([30;30;30]); %rad/sec
disp("Jacobian here:");
out2 = double(subs(J))
disp("Instantaneous joint velocities in mm/sec & rad/sec:");
inst = out2*joint_vel1

%5.c
F = [30*1000;0;0;0;0;0];%Newton= kgm/sec^2 *1000 =>mm
disp("Joint torques in Newton-mm:");
tau = out2'*F

%% Question 6 Inverse Velocity Kinematics

x_dot = [0;0;-100]; %mm/sec
Jv_inv = pinv(out2(1:3,:));
disp("Joint velocities in rad/sec:");
joint_vel2 = Jv_inv * x_dot

%% Question 7 Constraint Equation
syms alh beta t psi1 psi2 wc r a C xdot ydot tdot

%For alpha and beta values see figures above

%Rotation matrix from world to robot frame 
%t is theta = orientation of mobile robot
R = [cos(t), sin(t), 0;-sin(t), cos(t), 0;0,0,1]; %Rotmat I->R
eta_I = [xdot;ydot;tdot]; % in world frame

%Fixed standard wheel
roll = [sin(alh+beta),-cos(alh+beta),-(a/2)*cos(beta)];%rolling constraint
slide = [cos(alh+beta),sin(alh+beta),(a/2)*sin(beta)];%sliding constraint

%For left wheel
disp("Left wheel contraints");
%Rolling constraint 
%here 2*l = a
%wl = r*psi_dot / a
%wr = -r*psi_dot / a
%Therefore, r*psi_dot = wl*a (which is used below)
rc_left = roll * R* eta_I - r*psi1 
%Sliding constraint
sc_left = slide * R* eta_I

%Substituting alh, beta
alh = pi/2;beta = 0;
roll1 = (subs(roll));
slide1 = (subs(slide));
rc_left = roll1 * R* eta_I - r*psi1
sc_left = slide1 * R* eta_I


%For right wheel
disp("Right wheel constraints");
%Rolling constraint
rc_right = roll * R* eta_I + r*psi2
%Sliding constraint
sc_right = slide * R* eta_I

%Substituting alh, beta
alh = -pi/2; beta = pi;
roll2 = (subs(roll));
slide2 = (subs(slide));
rc_right = roll2 * R* eta_I + r*psi2
sc_right = slide2 * R* eta_I


%For castor wheel
rollc = [sin(alh+beta),-cos(alh+beta),(C)*cos(beta)];%rolling constraint
slidec = [cos(alh+beta),sin(alh+beta),C*sin(beta)];

disp("Castor wheel constraints");
%xr is given by = (r*psi_dot1 + r*psi_dot2)/2 

%Rolling constraint
rc_cast = rollc * R* eta_I - r*wc
%No Sliding constraint for omni-directional
sc_cast = slidec * R* eta_I 

%Substituting alh, beta
alh = pi; beta = pi/2;
roll3 = (subs(rollc));
slide3 = subs(slidec);
rc_cast = roll3 * R* eta_I - wc*r 
sc_cast = slide3 * R* eta_I

%% Question 8 Constraints Matrix
% eliminating redundant constraints of castor wheel
J_roll = [roll1;roll2;];
C_slide = [slide1];

%Combined constraints
J2 = [ r*psi1; -r*psi2];
zero = [0];
disp("Constraint matrix");
constraint = [J_roll;C_slide]
disp("Constraint equation");
eqn = constraint * R* eta_I == [J2;zero]



%% Question 9 Mobile Kinematics
syms xb yb %x,y of robot 
disp("Transformation Matrix from World to Reference Frame");
Tir = [R,[xb;yb;t];0,0,0,1]

%% Question 10 Numeric Solution
a = 507; r=143; psil= 30*2*pi*60; psi2=60*2*pi*60; C=300; %rad/sec
xb = 2.5*1000; yb= 1.5*1000; t = deg2rad(30);
out3 = double(subs(R));
disp("Instantaneous Velocities mm/sec");
eta_I = inv(out3) * [(psi1+psi2)/2;0;(psi1-psi2)/a]

%% Question 11 Combined Position Kinematics
T50 = inv(T04); %tip to robot
Tri = inv(Tir); %robot to world
disp("Transformation from tip to world/home frame");
MegaT = Tir * T04

%% Question 12 Numeric Solution
disp("Numeric Transformation:");
out4 = double(subs(MegaT))


%% Question 13 Combined Velocity Kinematics
disp("Tip velocities in world frame");
%Symbolic 6DOF linear and angular tip velocities
% syms lx ly lz ax ay az
% linear_world = MegaT(1:3,1:3) * [lx;ly;lz]
% angular_world = MegaT(1:3,1:3) * [ax;ay;az]

vel_world_base = Tir(1:3,1:3) * Jv*[t1;t2;t3];
total_vel = vel_world_base +eta_I;

xt_combined = total_vel;
disp("Combined Jacobian:");
syms psi1 psi2 t1 t2 t3   
J_combined = jacobian(xt_combined,[t1,t2,t3,psi1,psi2])

% Matrix vector equation
vel = J_combined * [dt1; dt2; dt3; psil; psi2];
%% Question 14 Force Propagation

b=424 ;c=300 ;D=380 ;e=328 ;f=323 ;g=82.4 ;
t1=deg2rad(20);t2=deg2rad(90);t3=deg2rad(30);
a = 507; r=143; psi1= 30*2*pi*60; psi2=60*2*pi*60; %rad/sec
xb = 2.5*1000; yb= 1.5*1000; t = deg2rad(30);

disp("Jacobian here:");
out5 = double(subs(J_combined))
disp("Torques required:");
tau_wheel = out5' * (out3(1:3,1:3)*[30*1000;0;0]) 
