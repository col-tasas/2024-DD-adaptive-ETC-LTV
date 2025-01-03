% generating A(k), B(k) as prescribed for the different cases

function [A, B] = sys_mats(CASE,k,T,mag_pert,T_delta,A_interpolated,B_interpolated)

% common LTI part for CASE=1,2
A_0 = [1.1 .1;
    0.1 0.2];

B_0 = [0.5 1;
    0.1 0.2];

switch CASE



    case 1


        A=A_0;

        B_1=B_0;

        B_2 = [B_1(1,1) -mag_pert*B_1(1,2);
            B_1(2,1) -mag_pert*B_1(2,2)];



        if ~mod(floor(abs(k-1)/T),2)



            B=B_1;

        else

            B=B_2;

        end



    case 2


        if T_delta~=0

            b=1; a=-1/T_delta;
            mag_pert_0=mag_pert;
            mag_pert = mag_pert_0*(a*k+b);

        end

        Pert=diag([mag_pert*cos(2*pi*k/T) -mag_pert*cos(2*pi*k/T)]);

        A = A_0*(eye(2)+Pert);

        B = B_0;



    case 3 % cubic plant from ODDAC paper (generated in file: AkBk_ODDAC)

A=A_interpolated(:,:,k);
 B=B_interpolated(:,:,k);



end



