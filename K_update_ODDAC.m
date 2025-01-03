function [K_0_ev, P_0_ev,Q_0_ev] = K_update_ODDAC(U_0,X_0,X_1,sigma_1,sigma_2,lambda,T_episode,LL)
% Algorithm proposed in
% S. Liu, K. Chen, and J. Eising, “Online data-driven adaptive control for
% unknown linear time-varying systems,” in IEEE Conference on Decision
% and Control, pp. 8775–8780, 2023.




% Extract data
T=size(X_0,2); % length of window where data are collected
nx=size(X_0,1); nu=size(U_0,1);

% Initialize SDP

Const = [];

alpha_1=sdpvar(1);alpha_2=sdpvar(1);


Q_i=sdpvar(nx);
L_i=sdpvar(nu,nx);

% Constraints and LMIs
Const =[ Const, alpha_1 >= 0]; 

Const =[ Const, alpha_2 >= 0]; 

Const =[ Const, Q_i >= (sigma_2^-1).* eye(nx)]; 
Const =[ Const, -Q_i >= -(sigma_1^-1).* eye(nx)]; % 10^(-3)




bar_M= [lambda.*Q_i zeros(nx,3*nx+2*nu);
zeros(nx,nx) zeros(nx,2*nx+2*nu) Q_i;
zeros(nu,nx) zeros(nu,2*nx+2*nu) L_i;
zeros(nx,nx) zeros(nx,2*nx+2*nu) Q_i;
zeros(nu,nx) zeros(nu,2*nx+2*nu) L_i;
zeros(nx,nx) Q_i L_i' Q_i L_i' Q_i]; 

pi_coeff=0;


for kkk=1:T

pi_coeff = pi_coeff + LL^2 * kkk^2 *norm( [X_0(:,kkk); U_0(:,kkk)] ,2)^2;

end

PI=pi_coeff.*eye(nx);

bar_N1_temp = [eye(nx) X_1;
          zeros(nx,nx) -X_0;
          zeros(nu,nx) -U_0;
          zeros(2*nx+nu,nx+T)];

bar_N1 = bar_N1_temp * [PI zeros(nx,T); zeros(T,nx) -eye(T)] * bar_N1_temp';


bar_N2_temp = [eye(nx) zeros(nx) zeros(nx,nu);
          zeros(nx+nu,2*nx+nu);
          zeros(nx) eye(nx) zeros(nx,nu);
          zeros(nu,2*nx) eye(nu);
          zeros(nx,2*nx+nu)];

bar_N2 = bar_N2_temp * [LL^2*T_episode^2*eye(nx) zeros(nx,nx+nu);
                        zeros(nx) -eye(nx) zeros(nx,nu); zeros(nu,2*nx) -eye(nu)] * bar_N2_temp';


LMI2 = bar_M-alpha_1*bar_N1-alpha_2*bar_N2;

Const =[ Const, LMI2 >= 0];





%obj=0; 
obj=trace(Q_i); 

options = sdpsettings('solver','mosek','verbose',0);

diagnostics = optimize(Const,obj,options);  %

if diagnostics.problem == 0 %%% feasible



    % fprintf('\n LMI design is feasible \n');



    Q_0val=value(Q_i);P_0val=inv(Q_0val);
        L_0val=value(L_i);
        K_0val=L_0val*P_0val;
        alpha_2_val=value(alpha_2);

    

K_0_ev=K_0val; 
P_0_ev=P_0val; 
Q_0_ev=Q_0val; 

else

    % fprintf('\n LMI design is infeasible \n');


    K_0_ev=[]; P_0_ev=[];  Q_0_ev=[];

end




end


