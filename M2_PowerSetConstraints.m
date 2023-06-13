function [demand,generation,fval,SelectionTime,OptimizationTime] = M2_PowerSetConstraints(T,N,A,P,E,Di,Gmax,Gmin)
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
% 2^T - 1 sub time intervals
W = [];
for t = 1:T
    Ind = nchoosek(1:T,t);
    WW = zeros(size(Ind,1),T);
    for tt = 1:size(Ind,1)
        WW(tt,Ind(tt,:)) = 1;
    end
    W = [W; WW];
end

% declare variables
d = sdpvar(T,1);
G = sdpvar(T,1);

% objective function
Objective = G'*G;

% constraints
Constraints = d >= 0;

% total energy
Constraints = [Constraints, ones(1,T)*d == ones(1,N)*E]; 

% real time constraints
for t = 1:T
    % power balance
    Constraints = [Constraints, G(t) == Di(t)+d(t)];
    
    % generator constraints
    Constraints = [Constraints, G(t) >= Gmin, G(t) <= Gmax];
end

% necessary and sufficient conditions
for i = 1:size(W,1)
    RHS = 0;
    for j = 1:N
        RHS = RHS + min(sum(A(j,:).*W(i,:)),E(j)/P(j))*P(j);
    end
    
    Constraints = [Constraints, W(i,:)*d <= RHS];
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
generation = value(G);
fval = value(Objective);

end

