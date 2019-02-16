clc
clear all

rng(345552)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                Chong CAO @ UCLA                  %%%
%%%   MGMTMFE405 Computational Methods in Finance    %%%
%%%            Professor Levon Goukasian             %%%
%%%             basis function (1,x,x^2)             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% American Option
r=0.06;
S0=38;
K=40;  % strike price
sig=0.40;
T=1.0;
N=100;   % number of time steps
Nsample=100000;   % number of samples


[Put,CI,Index] = AmericanOption(S0,K,r,T,sig,@PayOff,N,Nsample);
Put
CI


function [Price,PriceCI,Index] = AmericanOption(S0,K,r,T,sig,Fun,N,Nsim) 

dt=T/N;
Paths=pricepaths(S0,r,T,sig,N,Nsim);
EV=Fun(Paths,K);       % Exercise Value == PayOff
Index=zeros(Nsim,N+1); % matrix Index: 1 to exercise option, 0 to continue
                       % Each row can have at most one 1

% at time t=T
I=find(EV(:,N+1));  % find nonzero payoff
Index(I,N+1)=1;


% at time t=dt to t=T-dt
for i=1:N-1
    n=N+1-i;
    k=find(EV(:,n));   % find nonzero payoff
    X=EV(k,n);
    S=Paths(k,n);   % stock price
    
    M=length(X);
    Y=zeros(M,1);   % Continuation Value
    for j=1:M
        CV=0;
        for jj=n+1:N+1
            CV=CV+Index(k(j),jj)*EV(k(j),jj)*exp(-r*dt*(jj-n));
        end
    Y(j)=CV;
    end
    
    p=polyfit(S,Y,2);   % polynomial coefficients
    ECV=polyval(p,S);   % Expected Continuation Value
    
    % compare EV and ECV   
    for m=1:M
        if X(m)>ECV(m)
            Index(k(m),n)=1;
            Index(k(m),n+1:N+1)=0;
        else
            Index(k(m),n)=0;
        end
    end
    
end


Value=zeros(Nsim,0);    % option value == discounted payoff
for i=1:Nsim
    V=0;
    for j=2:N+1
        V=V+Index(i,j)*EV(i,j)*exp(-r*dt*(j-1));
    end
Value(i)=V;
end

% Mean value and Standard deviation
% To construct 95% confidence interval
[muhat,sigmahat,muci,sigmaci]=normfit(Value);
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
function Y = PayOff(X,K)

CallOrPut=0;   % 1 for call, 0 for put

if CallOrPut
    Y = max(X-K,0);
else
    Y = max(K-X,0);
end

end


