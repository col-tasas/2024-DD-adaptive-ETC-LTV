% hybrid systems data-based adaptive control core routine
function HS_DAC(CASE,CONTROLLER)

% Define parameters of time-varying perturbations

switch CASE

    case 1
        %% Switching plant

        % Specify:

        T_online = 100; % total simulated time

        % - period of the switch
        period= 12; % option analyzed in the paper: 12

        % - amplitude of the perturbation (input matrix second column)
        magnitude = 1; % options analyzed in the paper: 1; 2.5


        % options needed for the comparison with time-triggered controller (CONTROLLER=2)

        % period of the update
        trigger_period=8; % options analyzed in the paper: 8;12;16 (for mag=1)

        % offset in the triggering instant
        bias_T=0; % option analyzed in the paper: 0 (\neq 0 updates the controller with an offset compared to the switch)


        % not required but defined for compatibility of call to "sys_mats"
        T_delta = 0;
        A_interpolated=[]; B_interpolated=[];

    case 2
        %% Sinusoidal variations

        % Specify:

        T_online = 100; % total simulated time

        % - period of the sinusoid
        period= 10; % options analyzed in the paper: 8,10,12

        % - amplitude of the perturbation
        magnitude = 0.8; % option analyzed in the paper: 0.8
        T_delta = 30; % set to value \neq 0 if perturbation should decay to zero in finite time T_delta (last case in Figure 2)


        % (comparison with time-triggered controller  (CONTROLLER=2) not presented in the paper for this CASE, but comparison is possible in principle)

        % period of the update
        trigger_period=15;

        % offset in the triggering instant
        bias_T=0;


        % not required but defined for compatibility of call to "sys_mats"
        % T_delta = 0;
        A_interpolated=[]; B_interpolated=[];

    case 3
        %% Cubic plant from ODDAC paper
        % S. Liu, K. Chen, and J. Eising, “Online data-driven adaptive control for
        % unknown linear time-varying systems,” in IEEE Conference on Decision
        % and Control, pp. 8775–8780, 2023.

        [A_interpolated, B_interpolated] = AkBk_ODDAC;

        % Specify:
        T_online = 1000; % total simulated time

        scale_L=1.; % all entries of matrix A are scaled for all k by this value to test robustness of the controllers
        % options analyzed in the paper: 1, 1.1, 1.15, 1.2

        A_interpolated = scale_L.*A_interpolated;




        % (comparison with time-triggered controller  (CONTROLLER=2) not analyzed in the paper, but comparison is possible in principle)

        % period of the update
        trigger_period=15; %

        % offset in the triggering instant
        bias_T=0; %

        % not required but defined for compatibility of call to "sys_mats"
        T_delta = 0;
        period= 10; %
        magnitude = 0.8; %


end

% Getting nx and nu to initialize variables
k = 1; %dummy choice, just to get matrices sizes
[A_0, B_0] = sys_mats(CASE,k,period,magnitude,T_delta,A_interpolated,B_interpolated);
nx=size(A_0,1); nu=size(B_0,2);


% Offline data collection parameters to initialize data matrices
T_off = nx + nu ; % length of offline design (same as #columns of data-matrices online)

columns_data = T_off ;


if CASE==3
    x_0=ones(nx,1); % (initial condition as used in the ODDAC paper)
else

    x_0=zeros(nx,1); % zero initial condition
end

%% Offline experiment (starting at k=0)

if CONTROLLER==3 && CASE==3
    % ODDAC algorithm (only used for CASE==3) does not run offline phase and
    % uses pre-computed stabilizing K_0

    T_p = 30; % period of updates (ie the controller is updated every T_p timesteps): the value is taken from the shared working ODDAC implementation, it differs from the one reported in the ODDAC CDC paper which is 100

    columns_data = 10; % columns of data-matrices (from ODDAC paper specification)



    % Values taken from ODDAC paper
    K_0 =[0.13 0.26 -0.25 0.04 -0.13;
        0.08 0.28 0.13 0.05 0.01];


    x1=x_0; % (initial condition of the online phase, as used in the ODDAC paper)
    T_off = 0;


    % parameters ODDAC algorithm from CDC paper

    sigma_1=0.001; sigma_2=1000; lambda=0.9;

    LL = 0.003753345926796; % Lipschitz constant of time-varying matrices used for the robust design of K

    itt=1; % counter of time-triggered episode

    dith_v=10^-10; % magnitude of dithering signal used to excite system (in addition to feedback gain)


    % Initialize data-matrices

    X_0_ON=zeros(nx,columns_data); X_1_ON=X_0_ON; U_0_ON=zeros(nu,columns_data);



else

        m_off = 1; % magnitude of the offline input
    u_off = m_off .* rand(nu,T_off); % random uniform excitation


    % Initialize data-matrices
    X_0=zeros(nx,columns_data); X_1=X_0; U_0=zeros(nu,columns_data);

    x = x_0;

    for k = 1 : T_off

        u = u_off(:,k);
        x_h(:,k)=x; u_h(:,k)=u;

        [A, B] = sys_mats(CASE,k,period,magnitude,T_delta,A_interpolated,B_interpolated);

        A_h{k}=A;
        B_h{k}=B;

        x1 = A * x + B * u;

        X_0(:,k)=[x];
        X_1(:,k)=[x1];
        U_0(:,k)=[u];

        x=x1;

    end


    % First gain designed based on the open-loop collected data
    [K_0, S_0, F_0, a_1_0, a_2_0] = Map_L(U_0,X_0,X_1);


    % define online data-matrices (initialized at the off-line value)
    X_0_ON(:,1:columns_data)=X_0; X_1_ON(:,1:columns_data)=X_1; U_0_ON(:,1:columns_data)=U_0;


    % parameters event-triggered algorithm
    cnt_trig=0;
    simga_tol=0.1;

    S = S_0; a=a_1_0;
    sigma=1-simga_tol*(1-a); % needed for triggering condition
    tol_state = 3*10^(-3); % to avoid numerical issues, triggering only takes place when state norm is greater than this value


end

%% Online phase

% triggering condition

switch CONTROLLER

    case 0

        disp('Non-adaptive offline controller')

    case 1

        disp('Event-triggered controller')

    case 2

        disp('Time-triggered controller')

    case 3

        disp('ODDAC controller')

end



x=x1; % initial condition online phase
K = K_0;


for k = T_off + 1 : T_online % finite-horizon where problem is solved (starting from the last offline instant until T_online)



    if CONTROLLER==3

        i_new=floor(k/T_p);

        if k-i_new*T_p >=T_p-columns_data

            u = K*x + dith_v.*rand(nu,1); % dithering signal applied before planned controller update

            x_h(:,k)=x; u_h(:,k)=u;


            [A, B] = sys_mats(CASE,k,period,magnitude,T_delta,A_interpolated,B_interpolated);
            x1 = A * x + B * u;

            % saving data matrices

            X_0_ON_old=X_0_ON;X_1_ON_old=X_1_ON;U_0_ON_old=U_0_ON;
            X_0_ON(:,1)=[x]; X_1_ON(:,1)=[x1]; U_0_ON(:,1)=[u];
            X_0_ON(:,2:end)=X_0_ON_old(:,1:end-1); X_1_ON(:,2:end)=X_1_ON_old(:,1:end-1); U_0_ON(:,2:end)=U_0_ON_old(:,1:end-1);


        else

            u = K*x;
            x_h(:,k)=x; u_h(:,k)=u;
            [A, B] = sys_mats(CASE,k,period,magnitude,T_delta,A_interpolated,B_interpolated);
            x1 = A * x + B * u;


        end

    else

        u = K*x;

        x_h(:,k)=x; u_h(:,k)=u;

        [A, B] = sys_mats(CASE,k,period,magnitude,T_delta,A_interpolated,B_interpolated);

        A_h{k}=A; B_h{k}=B;


        x1 = A * x + B * u;

        % saving old data matrices
        X_0_ON_old=X_0_ON;X_1_ON_old=X_1_ON;U_0_ON_old=U_0_ON;

        % updating data matrices with new data
        X_0_ON(:,1)=[x]; X_1_ON(:,1)=[x1]; U_0_ON(:,1)=[u];
        X_0_ON(:,2:columns_data)=X_0_ON_old(:,1:end-1); X_1_ON(:,2:columns_data)=X_1_ON_old(:,1:end-1); U_0_ON(:,2:columns_data)=U_0_ON_old(:,1:end-1);

    end

    % Controller update


    if CONTROLLER==1

        % Checking event (Lyapunov condition not satisfied and state norm above threshold)

        if x1'*S*x1-sigma*x'*S*x >0 && norm(x1,2)>tol_state



            cnt_trig = cnt_trig + 1;
            idx_trig(cnt_trig,1)=k;


            [K_i, S_i, F_i, a_1_i, a_2_i] = Map_L(U_0_ON,X_0_ON,X_1_ON);

            if ~isempty(K_i)
                feas_design(cnt_trig,1)=1;

                K=K_i;S=S_i;a=a_1_i; sigma=1-simga_tol*(1-a_1_i); F=F_i;
                a_hist(cnt_trig,1)=a;
                S_norm(cnt_trig,1)=norm(S,2);


            else % map L is empty. Controller is not updated

                feas_design(cnt_trig,1)=0;

            end
            K_h{cnt_trig}=K;


        end


    elseif CONTROLLER==2

        % Periodically-triggered update

        % A new controller is designed with period "trigger_period" and bias
        % "bias_T" which represents an offset in controller update compared to when the plant actually switched


        if k>trigger_period && (floor( (k-(bias_T+columns_data+1))/trigger_period) >= (k-(bias_T+columns_data+1))/trigger_period) %&& norm(x1,2)>tol_state %&& norm(x1,2)>tol_state


            cnt_trig = cnt_trig + 1;
            idx_trig(cnt_trig,1)=k;
            [K_i, S_i, F_i, a_1_i, a_2_i] = Map_L(U_0_ON,X_0_ON,X_1_ON);

            if ~isempty(K_i)
                K=K_i;S=S_i;a=a_1_i;
                feas_design(cnt_trig,1)=1;
            else
                feas_design(cnt_trig,1)=0;

            end
            K_h{cnt_trig}=K;

        end




    elseif CONTROLLER==3

        i_new=floor(k/T_p);


        if itt==i_new % time to update controller

            scale = 1/norm(U_0_ON);     % scaling data-matrices for better numerical conditioning (from ODDAC implementation)

            X_0_ON_s=scale.*X_0_ON;
            X_1_ON_s=scale.*X_1_ON;
            U_0_ON_s=scale.*U_0_ON;


            [K_i, ~, ~] = K_update_ODDAC(U_0_ON_s,X_0_ON_s,X_1_ON_s,sigma_1,sigma_2,lambda,T_p,LL);


            if ~isempty(K_i)

                K=K_i;
                feas_design(itt,1)=1;
                K_h{itt}=K;


            else

                feas_design(itt,1)=0;

            end


            X_0_ON=zeros(nx,columns_data); X_1_ON=X_0_ON; U_0_ON=zeros(nu,columns_data);


            itt=itt+1;

        end



    end




    x=x1;

end



% quick plot showing first state component

time=1:k;
figure(1)
plot(time,x_h(1,time),'k'); hold on;

ylabel({'$x_1$'},'Interpreter','Latex')
xlabel({'$k$'},'Interpreter','Latex')


end





