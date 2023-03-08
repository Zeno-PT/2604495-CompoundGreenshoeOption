function P = AmericanTrinomial(S0, r, sigma, K, I, T, N, optiontype)

% This function compute the American Green shoe call/put
% option under the trinomial tree using a vector version.
% 
% INPUT:  S0 = initial stock price
%         r = annualized continuously compounded interest rate
%         sigma = annualized volatility of stock return
%         K = strike price
%         I = net price
%         T = time to expiration date of the option
%         N = number of time steps in the CRR model
%         optiontype = 'C' for call option, 'P' for put option
%
% OUTPUT: P = option price at time 0
%
% USAGE: P = AmericanTrinomial(40, 0.04, 0.3, 40, 39.5, 1, 10000, 'C')

dt = T/N;
u = exp(sigma*sqrt(2*dt));
% m = 1;
% d = 1/u;

u1 = exp(sigma*sqrt(dt/2));
d1 = 1/u1;

% calculate risk-neutral probability (pu, pm, pd)
pu = ((exp(r*dt/2)-d1)/(u1-d1))^2;
pd = ((u1-exp(r*dt/2))/(u1-d1))^2;
pm = 1-pu-pd;

% calculate option payoff
if strcmp(optiontype, 'C')
    Ovec = max(K-S0*u.^(N:-1:-N)', K-I); % 2N+1 nodes
elseif strcmp(optiontype,'P')
    Ovec = max(S0*u.^(N:-1:-N)'-K, K-I); % 2N+1 nodes
else
    error('"optiontype" has to be either ''C'' or ''P''.')
end

% calculate option price
for i=1:N
    n = (2*N+1)-(2*i); % n = number of nodes to be updated
    ContVal = exp(-r*dt)*(pu*Ovec(1:n, 1)+pm*Ovec(2:n+1, 1)+pd*Ovec(3:n+2, 1));
    if strcmp(optiontype, 'C')
        ExerVal = max(K-S0*u.^(N-i:-1:-N+i)', K-I);
    else
        ExerVal = max(S0*u.^(N-i:-1:-N+i)'-K, K-I);
    end
    if n==1 % no exercise is allowed at time 0
        ExerVal = 0;
    end
    Ovec(1:n, 1) = max(ExerVal, ContVal);
end

% get option price at time 0
P = Ovec(1,1);

end