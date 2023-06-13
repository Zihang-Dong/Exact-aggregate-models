function [demand,generation,fval,SelectionTime,OptimizationTime] = M3_ReducedConstraints_CombinatorialSearch(T,N,A,P,E,Di,Gmax,Gmin)
%%

% optimization
% options = sdpsettings('solver','fmincon','verbose',0);
% options.fmincon.Algorithm = 'sqp';
% options.fmincon.MaxIter = 10^(5);
% options.fmincon.MaxFunEvals = 10^(6);
% options.fmincon.TolCon = 1.0000e-12;
% options.fmincon.TolFun = 1.0000e-12;
% options.fmincon.TolFunValue = 1.0000e-12;

options = sdpsettings('solver','gurobi','verbose',0);

%% Iterative procedure of active constraint selection

WK = [];
Wk = ones(1,T);
AveCost_0 = Inf;
AveCost_History = [];

SelectionStart = tic;

for t = T:-1:1
    
    TSet = nchoosek(1:T,t);
    N_row = size(TSet,1);
    WW = zeros(N_row,T);
    
    if N_row > T
        Row = 1:1:N_row;
        for Column = 1:t
            WW(sub2ind([N_row,T],Row,TSet(:,Column)')) = 1;
        end
    else
        for tt = 1:N_row
            WW(tt,TSet(tt,:)) = 1;
        end
    end

    D_Sum = WW*Di';
    
    if N_row > 100000
        for j = 1:N
            D_Sum = D_Sum + min([WW*A(j,:)' E(j)/P(j)*ones(size(WW,1),1)],[],2)*P(j);
        end
    else
        N0 = min(N,1000);
        for k = 1:N/N0
            D_Sum_Temp = min(WW*A((k-1)*N0+1:(k-1)*N0+N0,:)'.*kron(ones(N_row,1),P((k-1)*N0+1:(k-1)*N0+N0)'), kron(ones(N_row,1),E((k-1)*N0+1:(k-1)*N0+N0)'));
            D_Sum = D_Sum + sum(D_Sum_Temp,2);
        end
    end
    
    
    AveCost = D_Sum./sum(WW,2);
    
    Wk = WW(AveCost == min(AveCost),:);
    
    if min(AveCost) < AveCost_0
        AveCost_0 = min(AveCost);
        AveCost_History = [AveCost_History; AveCost_0];
        WK = [WK; Wk];
    end
end

SelectionTime = toc(SelectionStart);

%% reduced constraint optimization problem

% declare variables
dd = sdpvar(T,1);
GG = sdpvar(T,1);

% objective function
Objective = GG'*GG;

% constraints
Constraints = dd >= 0;

% total energy
Constraints = [Constraints, ones(1,T)*dd == ones(1,N)*E]; 

% real time constraints
for t = 1:T
    % power balance
    Constraints = [Constraints, GG(t) == Di(t)+dd(t)];
    
    % generator constraints
    Constraints = [Constraints, GG(t) >= Gmin, GG(t) <= Gmax];
end

% necessary and sufficient conditions
rhs = zeros(size(WK,1),1);
for i = 1:size(WK,1)
    
    for j = 1:N
        rhs(i) = rhs(i) + min(sum(A(j,:).*WK(i,:)),E(j)/P(j))*P(j);
    end
    
    Constraints = [Constraints, WK(i,:)*dd <= rhs(i)];
end

% optimization
OptimizationStart= tic;
diagnostics = optimize(Constraints,Objective,options);
OptimizationTime = toc(OptimizationStart);

if diagnostics.problem ~= 0
    error('Something else happened')
end

%% optimal solutions
demand = value(dd);
generation = value(GG);
fval = value(Objective);

end

