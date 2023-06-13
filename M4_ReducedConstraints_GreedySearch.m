function [demand,generation,fval,SelectionTime,OptimizationTime] = M4_ReducedConstraints_GreedySearch(T,N,A,P,E,Di,Gmax,Gmin)
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

Wk = ones(1,T);
WK = Wk;
card = T;
WK_cost = Di*WK';
for j = 1:N
    WK_cost = WK_cost + min(sum(A(j,:).*WK),E(j)/P(j))*P(j);
end
WK_cost = WK_cost/T;
Wk_AveCost = Inf;

SelectionStart = tic;

while card > 1
    
    for i = card-1:-1:1

        Wk_Temp = [];
        
        for n = 1:size(Wk,1)
            Wkk = Wk(n,:);

            TSet = find(Wkk == 1);

            TSubSet = nchoosek(TSet,i);
            N_row = size(TSubSet,1);
            Wk0 = zeros(N_row,T);
            
            if N_row > T
                Row = 1:1:N_row;
                for Column = 1:i
                    Wk0(sub2ind([N_row,T],Row,TSubSet(:,Column)')) = 1;
                end
            else
                for tt = 1:N_row
                    Wk0(tt,TSubSet(tt,:)) = 1;
                end
            end
            
            D_Sum = Wk0*Di';
            
            if N_row > 100000
                for j = 1:N
                    D_Sum = D_Sum + min([Wk0*A(j,:)' E(j)/P(j)*ones(size(Wk0,1),1)],[],2)*P(j);
                end
            else
                N0 = min(N,1000);
                for k = 1:N/N0
                    D_Sum_Temp = min(Wk0*A((k-1)*N0+1:(k-1)*N0+N0,:)'.*kron(ones(N_row,1),P((k-1)*N0+1:(k-1)*N0+N0)'), kron(ones(N_row,1),E((k-1)*N0+1:(k-1)*N0+N0)'));
                    D_Sum = D_Sum + sum(D_Sum_Temp,2);
                end
            end

            Ave_Cost = D_Sum./sum(Wk0,2);

            if min(Ave_Cost) < Wk_AveCost
                Wk_AveCost = min(Ave_Cost);
                
                Wk_Temp = Wk0(Ave_Cost == min(Ave_Cost),:);
            elseif min(Ave_Cost) == Wk_AveCost
                Wk_Temp = [Wk_Temp; Wk0(Ave_Cost == min(Ave_Cost),:)];
            else
                continue;
            end

        end

        Wk_Temp = unique(Wk_Temp,'row');
        
        if Wk_AveCost < WK_cost(end)
            Wk = Wk_Temp;
            WK = [WK; Wk_Temp];
            WK_cost = [WK_cost; Wk_AveCost];
            break;
        end
    end
    
    card = i;
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

