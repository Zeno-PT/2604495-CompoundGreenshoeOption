function [P, SDP, CIP, OP] = SimCompoundOptionPrice(S0, r, sigma, K1, K2, I, T1, T2, optiontype, M)

% This function computes the compound Green Shoe option price based on the 
% Monte-Carlo simulation approximation with Antithetic sampling.
%
% INPUT:  S0 = initial stock price
%         r = annualized continuously compounded interest rate
%         sigma = annualized volatility of stock return
%         K1 = strike price for compound option
%         K2 = strike price for underlying option
%         I = net price
%         T1 = time to maturity (years) for compound option
%         T2 = time to maturity (years) for underlying option
%         optiontype =  'CoC' for call on call option,
%                       'CoP' for call on put option,
%                       'PoC' for put on call option,
%                       'PoP' for put on put option
%         M = number of simulations
% 
% OUTPUT: P = option price
%         SDP = SD of P
%         CIP = 95%-CI of P
%         OP = underlying option price at T1
%
% USAGE: [P, SDP, CIP, OP] = SimCompoundOptionPrice(40, 0.04, 0.3, 10, 40, 39.5, 0.5, 1, 'CoC', 1000)

N = 100;

% check optiontype input
if ~(strcmp(optiontype, 'CoC') || strcmp(optiontype, 'CoP') || ...
        strcmp(optiontype, 'PoC') || strcmp(optiontype, 'PoP'))
    error('"optiontype" has to be ''CoC'' or ''CoP'' or ''PoC'' or ''PoP''.')
end

% check T1 and T2 input
if T1 >= T2
    error('"T1" has to be less than "T2".')
end
if T1 < 0
    error('"T1" has to be greater than or equal to 0.')
end
if T2 < 0
    error('"T2" has to be greater than or equal to 0.')
end

% initial discounted payoff vectors
disc_payoff = nan(M, 1);
disc_payoff_hat = nan(M, 1);

% simulate stock price at time T1
epsilon = normrnd(0, 1, M, 1);
ST1 = S0*exp((r-0.5*sigma^2)*T1+sigma*sqrt(T1)*epsilon);
ST1_hat = S0*exp((r-0.5*sigma^2)*T1+sigma*sqrt(T1)*(-epsilon));

% calculate underlying option payoff
for i=1:M
    if mod(i, M/5)==0
        disp(['Calculating simulation ', num2str(i), ' from ', num2str(M), '.'])
    end
    disc_payoff(i) = AmericanTrinomial(ST1(i), r, sigma, K2, I, T2-T1, N, optiontype(3));
    disc_payoff_hat(i) = AmericanTrinomial(ST1_hat(i), r, sigma, K2, I, T2-T1, N, optiontype(3));
end

% calculate average of underlying option price at T1
OP = mean(0.5*(disc_payoff_hat+disc_payoff));

% calculate compound option payoff
if strcmp(optiontype(1), 'C')
    disc_payoff = exp(-r*T1)*max(disc_payoff-K1, 0);
    disc_payoff_hat = exp(-r*T1)*max(disc_payoff_hat-K1, 0);
else
    disc_payoff = exp(-r*T1)*max(K1-disc_payoff, 0);
    disc_payoff_hat = exp(-r*T1)*max(K1-disc_payoff_hat, 0);
end

% calculate compound option price
[P, S, CIP] = normfit((disc_payoff+disc_payoff_hat)/2);
SDP = S/sqrt(M);

end
