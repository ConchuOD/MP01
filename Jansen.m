function [ResultAngles, i_E, iter_count] = Jansen(Theta0, PreviousAngles, i_E)
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
% Max number of iterations, does not need to be high
max_count = 15;
% Assigning previous values as the next guess
ResultAngles = zeros(1,10);
ResultAngles(1) = PreviousAngles(1);
ResultAngles(2) = PreviousAngles(2);
ResultAngles(3) = PreviousAngles(3);
ResultAngles(4) = PreviousAngles(4);
ResultAngles(5) = PreviousAngles(5);
ResultAngles(6) = PreviousAngles(6);
ResultAngles(7) = PreviousAngles(7);
ResultAngles(8) = PreviousAngles(8);
ResultAngles(9) = PreviousAngles(9);
ResultAngles(10) = PreviousAngles(10);
%% This is a section for equation pair 1
Theta = [ResultAngles(1);ResultAngles(3)];
for N = 1:max_count
    % Constructing matrix [F1;F2]
    F1 = [a0*cos(Theta0)+ab*cos(Theta(1))+38-bc*cos(Theta(2));a0*sin(Theta0)+ab*sin(Theta(1))+7.8-bc*sin(Theta(2))]; 
    % Constructing matrix [dF1_T1,dF1_T2, dF2_T1,dF2_T2]
    dF1 = [-ab*sin(Theta(1)),bc*sin(Theta(2));ab*cos(Theta(1)),-bc*cos(Theta(2))];
    % Computing next iterations of Theta values
    Theta = Theta - (inv(dF1)*F1);  
    % Checking for convergence
    if (abs(F1(1))<Tol && abs(F1(2))<Tol)
        break
    end
end
% Warn user if no convergance
if (abs(a0*cos(Theta0)+ab*cos(Theta(1))+38-bc*cos(Theta(2))) > Tol || abs(a0*sin(Theta0)+ab*sin(Theta(1))+7.8-bc*sin(Theta(2))) > Tol ) % if difference is never greater than tolerance then warn user via console (no convergance)
    warning('Newton Raphson did not converge')
    i_E = i_E+1;
end

% Count number of iterations
iter_count(1)=N;
% Assign Theta values to result array
ResultAngles(1) = Theta(1);
ResultAngles(3) = Theta(2);
%% This is a section for equation pair 2
Theta = [ResultAngles(2);ResultAngles(4)];
for N = 1:max_count
    % Constructing matrix [F1;F2]
    F2 = [bd*cos(Theta(1))+bc*cos(ResultAngles(3))-cd*cos(Theta(2)); bd*sin(Theta(1))+bc*sin(ResultAngles(3))-cd*sin(Theta(2))];
    % Constructing matrix [dF1_T1,dF1_T2, dF2_T1,dF2_T2]
    dF2 = [-bd*sin(Theta(1)), cd*sin(Theta(2));bd*cos(Theta(1)), -cd*cos(Theta(2))];
    % Computing next iterations of Theta values
    Theta = Theta - (inv(dF2)*F2);  
    % Checking for convergence
    if (abs(F2(1))<Tol && abs(F2(2))<Tol)
        break
    end
end
% Warn user if no convergance
if (abs(bd*cos(Theta(1))+bc*cos(ResultAngles(3))-cd*cos(Theta(2))) > Tol || abs(bd*sin(Theta(1))+bc*sin(ResultAngles(3))-cd*sin(Theta(2))) > Tol ) % if difference is never greater than tolerance then warn user via console (no convergance)
    warning('Newton Raphson did not converge')
    i_E = i_E+1;
end
% Count number of iterations
iter_count(2)=N;
% Assign Theta values to result array
ResultAngles(2) = Theta(1);
ResultAngles(4) = Theta(2);
%% This is a section for equation pair 4
Theta = [ResultAngles(7);ResultAngles(8)];
for N = 1:max_count
    % Constructing matrix [F1;F2]
    F4 = [38+a0*cos(Theta0)-cf*cos(Theta(1))+af*cos(Theta(2)); 7.8+a0*sin(Theta0)-cf*sin(Theta(1))+af*sin(Theta(2))];
    % Constructing matrix [dF1_T1,dF1_T2, dF2_T1,dF2_T2]
    dF4 = [cf*sin(Theta(1)),-af*sin(Theta(2));-cf*cos(Theta(1)),af*cos(Theta(2))];
    % Computing next iterations of Theta values
    Theta = Theta - (inv(dF4)*F4);  
    % Checking for convergence
    if (abs(F4(1))<Tol && abs(F4(2))<Tol)
        break
    end
end
% Warn user if no convergance
if (abs(38+a0*cos(Theta0)-cf*cos(Theta(1))+af*cos(Theta(2))) > Tol || abs(7.8+a0*sin(Theta0)-cf*sin(Theta(1))+af*sin(Theta(2))) > Tol ) % if difference is never greater than tolerance then warn user via console (no convergance)
    warning('Newton Raphson did not converge')
    i_E = i_E+1;
end
% Count number of iterations
iter_count(3)=N;
% Assign Theta values to result array
ResultAngles(7) = Theta(1);
ResultAngles(8) = Theta(2);
%% This is a section for equation pair 3. This is done out of order as we would like Theta7 (ResultAngles 7) first, in order to get a better approximation.
Theta = [ResultAngles(5);ResultAngles(6)];
for N = 1:max_count
    % Constructing matrix [F1;F2]
    F3 = [cd*cos(ResultAngles(4))+de*cos(Theta(1))-cf*cos(ResultAngles(7))-ef*cos(Theta(2)); cd*sin(ResultAngles(4))+de*sin(Theta(1))-cf*sin(ResultAngles(7))-ef*sin(Theta(2))];
    % Constructing matrix [dF1_T1,dF1_T2, dF2_T1,dF2_T2]
    dF3 = [-de*sin(Theta(1)),ef*sin(Theta(2));de*cos(Theta(1)),-ef*cos(Theta(2))];
    % Computing next iterations of Theta values
    Theta = Theta - (inv(dF3)*F3);  
    % Checking for convergence
    if (abs(F3(1))<Tol && abs(F3(2))<Tol)
        break
    end
end
% Warn user if no convergance
if (abs(cd*cos(ResultAngles(4))+de*cos(Theta(1))-cf*cos(ResultAngles(7))-ef*cos(Theta(2))) > Tol || abs(cd*sin(ResultAngles(4))+de*sin(Theta(1))-cf*sin(ResultAngles(7))-ef*sin(Theta(2))) > Tol ) % if difference is never greater than tolerance then warn user via console (no convergance)
    warning('Newton Raphson did not converge')
    i_E = i_E+1;
end
% Count number of iterations
iter_count(4)=N;
% Assign Theta values to result array
ResultAngles(5) = Theta(1);
ResultAngles(6) = Theta(2);
%% This is a section for equation pair 5
Theta = [ResultAngles(9);ResultAngles(10)];
for N = 1:max_count
    % Constructing matrix [F1;F2]
    F5 = [cf*cos(ResultAngles(7))+gf*cos(Theta(1))-cd*cos(ResultAngles(4))-de*cos(ResultAngles(5))-eg*cos(Theta(2));cf*sin(ResultAngles(7))+gf*sin(Theta(1))-cd*sin(ResultAngles(4))-de*sin(ResultAngles(5))-eg*sin(Theta(2))];
    % Constructing matrix [dF1_T1,dF1_T2, dF2_T1,dF2_T2]
    dF5 = [-gf*sin(Theta(1)), eg*sin(Theta(2));gf*cos(Theta(1)), -eg*cos(Theta(2))];
    % Computing next iterations of Theta values
    Theta = Theta - (inv(dF5)*F5);  
    % Checking for convergence
    if (abs(F5(1))<Tol && abs(F5(2))<Tol)
        break
    end
end
% Warn user if no convergance
if (abs(cf*cos(ResultAngles(7))+gf*cos(Theta(1))-cd*cos(ResultAngles(4))-de*cos(ResultAngles(5))-eg*cos(Theta(2))) > Tol || abs(cf*sin(ResultAngles(7))+gf*sin(Theta(1))-cd*sin(ResultAngles(4))-de*sin(ResultAngles(5))-eg*sin(Theta(2))) > Tol ) % if difference is never greater than tolerance then warn user via console (no convergance)
    warning('Newton Raphson did not converge')
    i_E = i_E+1;
end
% Count number of iterations
iter_count(5)=N;
% Assign Theta values to result array
ResultAngles(9) = Theta(1);
ResultAngles(10) = Theta(2);
end
