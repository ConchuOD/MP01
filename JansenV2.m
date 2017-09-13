function [ResultAngles,i_E,iter_count] = JansenV2(Theta0, PreviousAngles, i_E)
% [ResultAngles] = Jansen(N, Theta0, PreviousAngles, i_E)
% Takes in crank angle Theta0 and a 1x10 matrix with the result angles for
% the last version. i_E is the error count.
% Setting up constants
Tol = 1E-9;
a0 = 15;
ab = 50;
bd = 55.8;
de = 39.4;
eg = 65.7;
gf = 49;
ef = 36.7;
af = 61.9;
cf = 39.3;
cd = 40.1;
bc = 41.5;
% Counting number of iterations
iter_count = (zeros(1,5));
% Max number of iterations
max_count = 15;
% Assigning first guesses for angles based on previous value
Theta = zeros(10,1);
Theta(1) = PreviousAngles(1);
Theta(2) = PreviousAngles(2);
Theta(3) = PreviousAngles(3);
Theta(4) = PreviousAngles(4);
Theta(5) = PreviousAngles(5);
Theta(6) = PreviousAngles(6);
Theta(7) = PreviousAngles(7);
Theta(8) = PreviousAngles(8);
Theta(9) = PreviousAngles(9);
Theta(10) = PreviousAngles(10);
for N = 1:max_count
    F1 = [ % equations
        a0*cos(Theta0)+ab*cos(Theta(1))+38-bc*cos(Theta(3));    % F1
        a0*sin(Theta0)+ab*sin(Theta(1))+7.8-bc*sin(Theta(3));   % F2
        bd*cos(Theta(2))+bc*cos(Theta(3))-cd*cos(Theta(4));     % F3
        bd*sin(Theta(2))+bc*sin(Theta(3))-cd*sin(Theta(4));     % F4
        38+a0*cos(Theta0)-cf*cos(Theta(7))+af*cos(Theta(8));    % F5
        7.8+a0*sin(Theta0)-cf*sin(Theta(7))+af*sin(Theta(8));   % F6
        cd*cos(Theta(4))+de*cos(Theta(5))-cf*cos(Theta(7))-ef*cos(Theta(6));    % F7
        cd*sin(Theta(4))+de*sin(Theta(5))-cf*sin(Theta(7))-ef*sin(Theta(6));    % F8
        cf*cos(Theta(7))+gf*cos(Theta(9))-cd*cos(Theta(4))-de*cos(Theta(5))-eg*cos(Theta(10));  % F9
        cf*sin(Theta(7))+gf*sin(Theta(9))-cd*sin(Theta(4))-de*sin(Theta(5))-eg*sin(Theta(10))]; % F10
    dF1 = [ % differential terms
        -ab*sin(Theta(1)),0,bc*sin(Theta(3)),0,0,0,0,0,0,0;                 % dF1_T1 dF1_T2  dF1_T3 dF1_T4 dF1_T5 dF1_T6 dF1_T7  dF1_T8 dF1_T9 dF1_T10 
        ab*cos(Theta(1)),0,-bc*cos(Theta(3)),0,0,0,0,0,0,0;                 % dF2_T1 dF2_T2  ...
        0,-bd*sin(Theta(2)),-bc*sin(Theta(3)),cd*sin(Theta(4)),0,0,0,0,0,0; % dF3_T1 ...
        0,bd*cos(Theta(2)),bc*cos(Theta(3)),-cd*cos(Theta(4)),0,0,0,0,0,0;  % dF4_T1 ...
        0,0,0,0,0,0,cf*sin(Theta(7)),-af*sin(Theta(8)),0,0;                 % dF5_T1 ...
        0,0,0,0,0,0,-cf*cos(Theta(7)),af*cos(Theta(8)),0,0;                 % dF6_T1 ...
        0,0,0,-cd*sin(Theta(4)),-de*sin(Theta(5)),ef*sin(Theta(6)),cf*sin(Theta(7)),0,0,0;  % dF7_T1 ...
        0,0,0,cd*sin(Theta(4)),de*cos(Theta(5)),-ef*cos(Theta(6)),-cf*cos(Theta(7)),0,0,0;  % dF8_T1 ...
        0,0,0,cd*sin(Theta(4)),de*sin(Theta(5)),0,-cf*sin(Theta(7)),0,-gf*sin(Theta(9)), eg*sin(Theta(10));     % dF9_T1 ...
        0,0,0,-cd*cos(Theta(4)),-de*cos(Theta(5)),0,cf*cos(Theta(7)),0,gf*cos(Theta(9)),-eg*cos(Theta(10))];    % dF10_T1 ...
    % Assigning next value of Angles
    Theta = Theta - (inv(dF1)*F1);
    % Checking for convergence
    if (abs(F1(1))<Tol && abs(F1(2))<Tol && abs(F1(3))<Tol && abs(F1(4))<Tol && abs(F1(5))<Tol && abs(F1(6))<Tol && abs(F1(7))<Tol && abs(F1(8))<Tol && abs(F1(9))<Tol && abs(F1(10))<Tol)
       break
    end
end
% Warn user if did not converge
if ( (abs(F1(1))>Tol) || (abs(F1(2))>Tol) || (abs(F1(3))>Tol) || (abs(F1(4))>Tol) || (abs(F1(5))>Tol) || (abs(F1(6))>Tol) || (abs(F1(7))>Tol) || (abs(F1(8))>Tol) || (abs(F1(9))>Tol) || (abs(F1(10))>Tol) ) % if difference is never greater than tolerance then warn user via console (no convergance)
   warning('Newton Raphson did not converge')
   i_E = i_E+1;
end
ResultAngles = Theta';
end
