% Map L returning variables that satisfy Property P_T(d)
function [K_0, S_0, F_0val, a_1val, a_2val] = Map_L(U_0,X_0,X_1)


% Extract dimensions
T=size(X_0,2); nx=size(X_0,1);

% Initialize SDP

Const = [];


alpha_0=sdpvar(1);
Q_0=sdpvar(T,nx);
H_0=sdpvar(nx);
P_0=sdpvar(nx);

Const =[ Const, alpha_0 >= 0];

Const =[ Const, H_0 >= 0];

Const =[ Const, P_0 == X_0*Q_0];

LMI1 = [P_0-alpha_0.*X_1*(X_1')-H_0   X_1*Q_0;
    (X_1*Q_0)' P_0];


Const =[ Const, LMI1 >= 0];



LMI2 = [eye(T)  Q_0;
    Q_0' P_0];

Const =[ Const, LMI2 >= 0];

obj=-geomean(H_0); % maximizing the geomean ~determinant of H_0




options = sdpsettings('solver','mosek','verbose',0);


diagnostics = optimize(Const,obj,options);



if diagnostics.problem == 0



    % fprintf('\n LMI design is feasible \n');


    % needed for the controller triggering
    P_0val=value(P_0);
    S_0=inv(P_0val);
    Q_0val=value(Q_0);
    K_0=U_0*Q_0val*S_0;

    % needed for triggering

    alpha_0val=value(alpha_0);
    H_0val=value(H_0);
    alpha_frac=alpha_0val/(alpha_0val+1);


    epsilon_W=0.1;

    F_0val=alpha_frac*H_0val-epsilon_W*alpha_frac*H_0val;


    W_0 = H_0val - 1/alpha_frac * F_0val;

    a_0=sdpvar(1);
    Const2=[];
    Const2 =[ Const2, a_0 >= 0];
    Const2 =[ Const2, W_0 >= a_0* P_0val];


    obj2=-a_0; % maximizing a_0


    diagnostics2 = optimize(Const2,obj2,options);

    a_0_tilde = value(a_0);

    a_1val = 1-a_0_tilde;
    a_2val = 1+1/alpha_0val;



    if a_1val>1
        % keyboard
        fprintf('\n Warning: Numerical problems in computing a_0: set by default to 0.01 \n');
a_0_tilde = 0.01;

    a_1val = 1-a_0_tilde;
    a_2val = 1+1/alpha_0val;
    end




else

    % fprintf('\n LMI design is infeasible \n');

    K_0=[]; S_0=[]; a_1val=[]; F_0val=[]; a_2val=[];

end




end


