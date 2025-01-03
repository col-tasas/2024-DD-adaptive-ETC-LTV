% This function generates the time-varying matrices A(k),B(k) proposed in
% the CDC paper where ODDAC method used here for comparison was proposed:

% S. Liu, K. Chen, and J. Eising, “Online data-driven adaptive control for
% unknown linear time-varying systems,” in IEEE Conference on Decision
% and Control, pp. 8775–8780, 2023.

function [A_interpolated, B_interpolated] = AkBk_ODDAC



A_0=[-0.5 -0.4 0.1 -0.8 -0.2;
-0.5 -0.1 0.2 0.7 0;
-0.4 -0.9 0.6 -0.3 0.4;
0.2 -0.3 -1.2 0 -0.1;
-0.6 0.8 -0.5 -0.1 -0.1];

B_0=[-1.4 2.2;
0.9 1.4;
2.7 0.5;
-0.7 1.5;
0.6 -1.9];

A_1=[-0.5 -0.7 0.3 -0.6 0;
0 0 0 0.8 0.4;
-0.7 -1.0 0.7 0.1 0.2;
-0.2 -0.2 -1.1 0.3 0.3;
-0.9 0.7 -0.9 0.5 0.4];

B_1=[-1.5 2.4;
0.9 1.3;
2.9 0.7;
-0.7 1.5;
0.4 -1.9];

A_2=[0 -0.6 -0.2 -0.7 0.5
0 0.1 0.4 1.1 0.7
-1.4 -0.9 0.5 0.5 0.5
-0.2 -0.2 -1.5 -0.3 0.5
-0.9 0.5 -0.6 0.7 0.5];

B_2=[-1.4 2.4
0.9 1.5
3.0 0.6
-0.8 1.5
0.5 -1.9];



T=1000;
A_matrix3D = cat(3, A_0, A_1, A_2);
B_matrix3D = cat(3, [B_0 zeros(5,1)], [B_1 zeros(5,1)], [B_2 zeros(5,1)]); % adding a fictitious column to B to make interp3 work with cubic approximation despites B has only 2 columns...

% Define the original grid points
[x_A, y_A, z_A] = meshgrid(1:size(A_matrix3D,2), 1:size(A_matrix3D,1), 1:size(A_matrix3D,3));
[x_B, y_B, z_B] = meshgrid(1:size(B_matrix3D,2), 1:size(B_matrix3D,1), 1:size(B_matrix3D,3));
[xq_A, yq_A, zq_A] = meshgrid(1:1:size(A_matrix3D,2), 1:1:size(A_matrix3D,1), linspace(1,size(A_matrix3D,3),T)); % linspace(x1,x2,n)
[xq_B, yq_B, zq_B] = meshgrid(1:1:size(B_matrix3D,2), 1:1:size(B_matrix3D,1), linspace(1,size(B_matrix3D,3),T)); % linspace(x1,x2,n)


% Perform cubic interpolation
A_interpolated = interp3(x_A, y_A, z_A, A_matrix3D, xq_A, yq_A, zq_A, 'cubic');
B_interpolated_0 = interp3(x_B, y_B, z_B, B_matrix3D, xq_B, yq_B, zq_B, 'cubic');
B_interpolated=B_interpolated_0(1:end,1:end-1,1:end);


%% Test: plotting the 1,1 entry of A and B
%x = 1:1:size(A_matrix3D,3);
%y = [A_0(1,1) A_1(1,1) A_2(1,1)];
%plot(x,y,'o',linspace(1,size(A_matrix3D,3),T),squeeze(A_interpolated(1,1,:)))

% x = 1:1:size(B_matrix3D,3);
% y = [B_0(1,1) B_1(1,1) B_2(1,1)];
% plot(x,y,'o',linspace(1,size(B_matrix3D,3),T),squeeze(B_interpolated(1,1,:)))
% 
% y = [B_0(3,1) B_1(3,1) B_2(3,1)];
% plot(x,y,'o',linspace(1,size(B_matrix3D,3),T),squeeze(B_interpolated(3,1,:)))
% 
% 
% keyboard

end


% tests with spline
% x = [0 500 1000];
% y = [A_0(1,1) A_1(1,1) A_2(1,1)];
% xx = 0:.25:1000;
% yy = spline(x,y,xx);
% plot(x,y,'o',xx,yy)