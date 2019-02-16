clc
clear all

randn('seed',0)

% American Option
r=0.1;
S0=200;
K=220;
sig=0.30;
T=1;
N=250;   % number of time steps
Nsample=100000;   % number of samples


[Put,CI] = AmericanOption(S0,K,r,T,sig,@PayoffFun,N,Nsample)


function [Price,PriceCI] = AmericanOption(S0,K,r,T,sig,Fun,N,Nsim) 

dt=T/N;
Paths=pricepaths(S0,r,T,sig,N,Nsim);
PayOff=Fun(Paths,K);
V=zeros(Nsim,N+1);
V(:,N+1)=PayOff(:,N+1);   % value of option

for i=1:N-1
    n=N+1-i;    
    X=PayOff(:,n);
    k=find(X);   % find nonzero payoff
    XX=X(k);
    S=Paths(k,n);   % stock price
    Y=V(k,n+1)*exp(-r*dt);
    p=polyfit(S,Y,2);
    YY=polyval(p,S);
    
    % compare XX and YY
    M=length(XX);
    for m=1:M
        if XX(m)>YY(m)
            V(k(m),n)=PayOff(k(m),n);
        else
            V(k(m),n)=V(k(m),n+1)*exp(-r*dt);
        end
    end
    
    kk=find(~X);   % find zero payoff
    V(kk,n)=V(kk,n+1)*exp(-r*dt);
end

% Mean value and Standard deviation
% To construct 95% confidence interval
[muhat,sigmahat,muci,sigmaci]=normfit(V(:,2)*exp(-r*dt));
Price=muhat;
PriceCI=muci;

end


% Geometric Brownian Motion
% dX(t)=r*X(t)*dt+sig*X(t)*dW(t)
function Paths = pricepaths(S0,r,T,sig,N,Nsim)

DeltaT = T/N; 
Drift = (r-sig^2/2)*DeltaT; 
Volatility = sig*sqrt(DeltaT)*randn(Nsim,N); 
Increments = Drift+Volatility;
LogPaths = cumsum([log(S0)*ones(Nsim,1),Increments],2); 
Paths = exp(LogPaths);

end



% payoff function
function Y = PayoffFun(X,K)

CallOrPut=0;   % 1 for call, 0 for put

if CallOrPut
    Y = max(X-K,0);
else
    Y = max(K-X,0);
end

end
