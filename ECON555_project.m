% Elijah Adrian %

% In vectors and matrices, along columns and rows,
% the ordering is Canada, US, Mexico, Rest of the World (RoW)

% I take a very function-based approach to solving the model and make heavy
% use of MATLAB matrix operations for compactness

% SECTION 1 %
% Step 1 - placeholder values
L = [1; 10; 5; 100];
T = [1; 2; 0.5; 0.5];
d = [1 2 2 2; 2 1 2 2; 2 2 1 2; 2 2 2 1];

% Step 2 - more placeholder values
theta = 4;
tau = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % default no tariffs

% Step 3 - multilateral resistance with w_0
w0 = [1; 1; 1; 1];
function Phi = Phi_maker(w,T,d,theta,tau)
    d_tilde = d .* (1 + tau);
    Phi = (T'*((w .* d_tilde').^(-theta)))';
end
Phi = Phi_maker(w0,T,d,theta,tau)

% Step 4 - bilateral trade share matrix with w_0
function pi = pi_maker(w,T,d,theta,tau)
    d_tilde = d .* (1 + tau);
    pi = ((T.*((w .* d_tilde').^(-theta))) ./ (Phi_maker(w,T,d,theta,tau))')';
end
pi = pi_maker(w0,T,d,theta,tau)

% Step 5 - excess expenditure function and vector
function EXE = model_iteration(w,T,L,d,theta,tau)
    EXE = w .* L - (pi_maker(w,T,d,theta,tau)')*(w .* L);
end
EXE = model_iteration(w0,T,L,d,theta,tau)

% Step 6 - non-normalized equilibrium wages
fcf = @(w)model_iteration(w,T,L,d,theta,tau);
[w_equ] = fsolve(fcf,w0);

% Step 7 - normalized equilibrium wages and equilibrium excess expenditure
w_star = w_equ/w_equ(1) % normalizing to w_c = 1
model_iteration(w_star,T,L,d,theta,tau) % very close to zero, as expected

% Step 8 - equilibrium with higher US technology
T = [1; 3; 0.5; 0.5];
fcf = @(w)model_iteration(w,T,L,d,theta,tau);
[w_equ2] = fsolve(fcf,w0);
w_star2 = w_equ2/w_equ2(1) % w_u increases from 0.86 to 0.94


% SECTION 2 %
%Step 1 - read trade data
trade = readtable("X_ni-2023.csv", 'VariableNamingRule','preserve');
X_ni = str2double(table2array(trade(:,[2,3,4,5])));

% Step 2 - read labour and GDP data, then solve for and normalize wages and
% solve for bilateral trade shares
WBI = readtable("projectecon555wdidata.csv", 'VariableNamingRule','preserve');
L_data = table2array(WBI([2,6,4,8],[5]))

GDP_data = table2array(WBI([1,5,3,7],[5]));
w_data_nn = GDP_data ./ L_data
w_data = w_data_nn/w_data_nn(1)

X_ni(isnan(X_ni)) = 0;
for i = 1:4
    X_ni(i,i) = GDP_data(i)-ones([1,4])*X_ni(i,:)';
end
pi_data = X_ni ./ GDP_data % Mexico is the most open, RoW is least open

% Step 3 - solve for technology
T_data = (w_data(1)./w_data).^(-theta);
for i = 1:4
    T_data(i) = T_data(i) * (((pi_data(i,i)/pi_data(i,1))/(pi_data(1,1)/pi_data(1,i)))^(1/2));
end
T_data

% Step 4 - solve for trade costs
d_data = zeros([4,4]);
for n = 1:4
    for i = 1:4
        d_data(n,i) = (((T_data(n)*pi_data(n,i))/(T_data(i)*pi_data(n,n)))^(-1/theta))*(w_data(n)/w_data(i));
    end
end
d_data


% SECTION 3 %
% new function to solve for excess expenditure with tariffs
function EXE_tariff = model_iteration_tariff(w,T,L,d,theta,tau)
    sum = diag(pi_maker(w,T,d,theta,tau) * tau');
    Gamma = (w .* L) .* (sum ./ (1- sum));
    EXE_tariff = (w .* L) + Gamma - (pi_maker(w,T,d,theta,tau)')*((w .* L) + Gamma);
end

% new function to solve for equilibrium empirical wage vector with tariffs
function tariff_wage = tariff_wage_generator(w,T,L,d,theta,tau)
    fcf2 = @(w)model_iteration_tariff(w,T,L,d,theta,tau);
    [wage] = fsolve(fcf2,w);
    tariff_wage = wage/wage(1);
end

% empirical value of EXE with surpluses/deficits
S = model_iteration_tariff(w_data,T_data,L_data,d_data,theta,tau)

% find baseline (no tariffs) equilibrium empirical wages and pi without surpluses/deficits
wage_baseline = tariff_wage_generator(w_data,T_data,L_data,d_data,theta,tau);
pi_baseline = pi_maker(wage_baseline,T_data,d_data,theta,tau)

% function for calculating welfare changes against baseline
function welfare_change = welfare_change_func(w,T,L,d,theta,tau)
    tau_baseline = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    wage_baseline = tariff_wage_generator(w,T,L,d,theta,tau_baseline);
    pi_baseline = pi_maker(wage_baseline,T,d,theta,tau_baseline);
    wage_hat = tariff_wage_generator(w,T,L,d,theta,tau);
    pi_hat = pi_maker(wage_hat,T,d,theta,tau);
    U_hat = diag(pi_hat ./ pi_baseline) .^ (-1/theta);
    welfare_change = 100 .* (U_hat - 1);
end

% Question 1
% simulation 1a
tau_1a = [0 0 0 0; 0.25 0 0 0; 0 0 0 0; 0 0 0 0];
welfare_1a = welfare_change_func(w_data,T_data,L_data,d_data,theta,tau_1a)

% simulation 1b
tau_1b = [0 0 0 0; 0 0 0.25 0; 0 0 0 0; 0 0 0 0];
welfare_1b = welfare_change_func(w_data,T_data,L_data,d_data,theta,tau_1b)

% simulation 1c
tau_1c = [0 0 0 0; 0.25 0 0.25 0; 0 0 0 0; 0 0 0 0];
welfare_1c = welfare_change_func(w_data,T_data,L_data,d_data,theta,tau_1c)

% simulation 1d
tau_1d = [0 0.25 0 0; 0.25 0 0.25 0; 0 0 0 0; 0 0 0 0];
welfare_1d = welfare_change_func(w_data,T_data,L_data,d_data,theta,tau_1d)

% simulation 1e
tau_1e = [0 0 0 0; 0.25 0 0.25 0; 0 0.25 0 0; 0 0 0 0];
welfare_1e = welfare_change_func(w_data,T_data,L_data,d_data,theta,tau_1e)

% simulation 1f
tau_1f = [0 0.25 0 0; 0.25 0 0.25 0; 0 0.25 0 0; 0 0 0 0];
welfare_1f = welfare_change_func(w_data,T_data,L_data,d_data,theta,tau_1f)

% Question 2
tau_2 = [0 0 0 0; 0.25 0 0 0; 0 0 0 0; 0 0 0 0];

% new welfare change function incorporating reduction in trade costs
function welfare_change2 = welfare_change_func2(w,T,L,d,theta,tau,A)
    tau_baseline = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    d_hat = d;
    d_hat(1,4) = A * d_hat(1,4);
    d_hat(4,1) = A * d_hat(4,1);
    wage_baseline = tariff_wage_generator(w,T,L,d,theta,tau_baseline);
    pi_baseline = pi_maker(wage_baseline,T,d,theta,tau_baseline);
    wage_hat = tariff_wage_generator(w,T,L,d_hat,theta,tau);
    pi_hat = pi_maker(wage_hat,T,d_hat,theta,tau);
    U_hat = diag(pi_hat ./ pi_baseline) .^ (-1/theta);
    welfare_change2 = 100 .* (U_hat - 1);
end

% simulation 2a
A_a = 0.975;
welfare_2a = welfare_change_func2(w_data,T_data,L_data,d_data,theta,tau_2,A_a)

% simulation 2b
A_b = 0.95;
welfare_2b = welfare_change_func2(w_data,T_data,L_data,d_data,theta,tau_2,A_b)

% simulation 2c
A_c = 0.925;
welfare_2c = welfare_change_func2(w_data,T_data,L_data,d_data,theta,tau_2,A_c)

% simulation 2d
A_d = 0.9;
welfare_2d = welfare_change_func2(w_data,T_data,L_data,d_data,theta,tau_2,A_d)

% END OF CODE