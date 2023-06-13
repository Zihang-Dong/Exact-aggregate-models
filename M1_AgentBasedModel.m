function [demand,generation,fval,SelectionTime,OptimizationTime] = M1_AgentBasedModel(T,N,A,P,E,Di,Gmax,Gmin)
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

SelectionTime = 0;

%% Optimization problem

% declare variables
d = sdpvar(1,T);
G = sdpvar(1,T);
u = sdpvar(N,T);

% objective function
Objective = G*G';

% constraints
Constraints = d >= 0;

% real time constraints
for t = 1:T
    % power balance
    Constraints = [Constraints, G(t) == Di(t)+d(t)];
    
    % generator constraints
    Constraints = [Constraints, G(t) >= Gmin; G(t) <= Gmax];
    
    % aggregate demand
    Constraints = [Constraints, d(t) == sum(u(:,t))];
    
    % power
    for j = 1:N
        Constraints = [Constraints, 0 <= u(j,t) <= P(j).*A(j,t)];
    end
    
end

% energy constraint
for j = 1:N
    Constraints = [Constraints, u(j,:)*ones(T,1) >= E(j)];
end
        
% optimization problem
OptimizationStart= tic;
diagnostics = optimize(Constraints,Objective,options);
OptimizationTime = toc(OptimizationStart);

if diagnostics.problem ~= 0
    error('Something else happened')
end


%% optimal solutions
demand = value(d);
demand = demand';
generation = value(G);
generation = generation';
fval = value(Objective);

end

