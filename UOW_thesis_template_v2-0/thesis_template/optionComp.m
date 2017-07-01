function optionComp

%% Financial parameters %%
K = 100;
r = 0.05;
q = 0.00;
T = 0.5;
sigma = 0.3;
S0 = 20:0.5:500;

%% Jump parameters %%
lambda = 0.5;
mu = -0.90;
delta = 0.45;
iter = 100;
alpha = exp(mu + 0.5*delta^2) - 1;
lambdaF = lambda*(1+alpha);

%% Useful functions %%
N = @(x) 0.5 + 0.5*erf(x/sqrt(2));
z1 = @(x,t,u) (log(x) + (r-q+0.5*sigma^2).*(u-t)) ./ (sigma*sqrt(u-t));
z2 = @(x,t,u) (log(x) + (r-q-0.5*sigma^2).*(u-t)) ./ (sigma*sqrt(u-t));
call = @(x,t) -K.*exp(-r.*(T-t)).*N(z2(x./K,t,T)) ...
        + x.*exp(-q.*(T-t)).*N(z1(x./K,t,T));

%% Computing the jump-diffusion European call option price %%
Cj = zeros(length(S0),1);
for i = 1:length(S0)
    sum = 0;
    for n = 0:iter
        rF = r - lambda*alpha + (n*log(1+alpha))/T;
        sigmaF = sqrt(sigma^2 + (n*delta^2)/T);
        d1 = (log(S0(i)./K) + (rF - q + 0.5*sigmaF^2)*T)/(sigmaF*sqrt(T));
        d2 = d1 - sigmaF*sqrt(T);
        JF = ((lambdaF*T)^n/factorial(n))*exp(-lambdaF*T);
        CBSM = S0(i).*exp(-q*T)*N(d1) - K*exp(-rF*T)*N(d2);
        sum = sum + JF*CBSM;
    end
    Cj(i) = sum;
end
C = call(S0(:),0);

%% Plotting the results %%
set(gca,'FontSize',18,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','bold')

figure (1);
plot(S0,C,'r',S0,Cj,'b--');
legend('Black-Scholes','Lognormal jumps',4);
axis([0 500 0 max(Cj)]);
title('Black-Scholes versus Merton''s model with lognormal jumps');
xlabel('Asset price (S_0)');
ylabel('Call option value (V_0(S_0,0))');
grid on;

figure (2);
plot(S0,Cj-C);
title('Difference between Merton`s model and Black-Scholes');
xlabel('Asset price (S_0)');
ylabel('Call option difference (v_M - v)');
grid on;
end