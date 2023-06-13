clear all;
close all;
clc;
yalmip('clear')
warning off

% Number of scenarios
NS = 2;

% number of devices
NVec = [10 100 1000 10000];
% NVec = [10 20 30 40];
% NVec = 10000;

% time window
TVec = [8 16 24];
% TVec = [4 6 8];
% TVec = 24;

% struct
Obj(1:length(NVec)*length(TVec)) = struct;

% Implementation
for n = 1:length(NVec)
    N = NVec(n);
    for t = 1:length(TVec)
        T = TVec(t);
        
        % devices' power and energy
        P0 = ones(N,NS);
        E0 = zeros(N,NS);
        
        % devices' availability window
        A0 = zeros(N,T,NS);
        for s = 1:NS
            for i = 1:N
                C = nchoosek(1:T,randi([1 T]));
                A0(i,C(randi([1 size(C,1)]),:),s) = 1;
                E0(i,s) = rand(1)*(sum(A0(i,:,s)))*P0(i,s);
            end
        end
        
%         load('A0E0.mat');
        
        % inflexible demand
        Di0 = randi([1 5*(N/10)],NS,T);
        
        % Generation
        Gmax = 15*(N/10);
        Gmin = 0;
        
        % data
        Obj((n-1)*length(TVec)+t).T = T;
        Obj((n-1)*length(TVec)+t).N = N;
        Obj((n-1)*length(TVec)+t).A0 = A0;
        Obj((n-1)*length(TVec)+t).P0 = P0;
        Obj((n-1)*length(TVec)+t).E0 = E0;
        Obj((n-1)*length(TVec)+t).Gmax = Gmax;
        Obj((n-1)*length(TVec)+t).Gmin = Gmin;
        Obj((n-1)*length(TVec)+t).Di0 = Di0;
        Obj((n-1)*length(TVec)+t).FailureNumber_M4 = 0;
        Obj((n-1)*length(TVec)+t).Success_Rate_M4 = 1;
        Obj((n-1)*length(TVec)+t).FailureNumber_M5 = 0;
        Obj((n-1)*length(TVec)+t).Success_Rate_M5 = 1;
        Obj((n-1)*length(TVec)+t).demand_M1 = zeros(T,NS);
        Obj((n-1)*length(TVec)+t).generation_M1 = zeros(T,NS);
        Obj((n-1)*length(TVec)+t).demand_M2 = zeros(T,NS);
        Obj((n-1)*length(TVec)+t).generation_M2 = zeros(T,NS);
        Obj((n-1)*length(TVec)+t).demand_M3 = zeros(T,NS);
        Obj((n-1)*length(TVec)+t).generation_M3 = zeros(T,NS);
        Obj((n-1)*length(TVec)+t).demand_M4 = zeros(T,NS);
        Obj((n-1)*length(TVec)+t).generation_M4 = zeros(T,NS);
        Obj((n-1)*length(TVec)+t).demand_M5 = zeros(T,NS);
        Obj((n-1)*length(TVec)+t).generation_M5 = zeros(T,NS);
        Obj((n-1)*length(TVec)+t).SelectionTime_M1 = zeros(1,NS);
        Obj((n-1)*length(TVec)+t).OptimizationTime_M1 = zeros(1,NS);
        Obj((n-1)*length(TVec)+t).SelectionTime_M2 = zeros(1,NS);
        Obj((n-1)*length(TVec)+t).OptimizationTime_M2 = zeros(1,NS);
        Obj((n-1)*length(TVec)+t).SelectionTime_M3 = zeros(1,NS);
        Obj((n-1)*length(TVec)+t).OptimizationTime_M3 = zeros(1,NS);
        Obj((n-1)*length(TVec)+t).SelectionTime_M4 = zeros(1,NS);
        Obj((n-1)*length(TVec)+t).OptimizationTime_M4 = zeros(1,NS);
        Obj((n-1)*length(TVec)+t).SelectionTime_M5 = zeros(1,NS);
        Obj((n-1)*length(TVec)+t).OptimizationTime_M5 = zeros(1,NS);

        % dvaluation
        for s = 1:NS
            
            n
            t
            s
            
            % Method 1: Detailed agent based model
            [demand_M1,generation_M1,SelectionTime_M1,OptimizationTime_M1] = M1_AgentBasedModel(T,N,A0(:,:,s),P0(:,s),E0(:,s),Di0(s,:),Gmax,Gmin);
            Obj((n-1)*length(TVec)+t).demand_M1(:,s) = demand_M1;
            Obj((n-1)*length(TVec)+t).generation_M1(:,s) = generation_M1;
            Obj((n-1)*length(TVec)+t).SelectionTime_M1(:,s) = SelectionTime_M1;
            Obj((n-1)*length(TVec)+t).OptimizationTime_M1(:,s) = OptimizationTime_M1;
            fprintf('\nM1 finished! \n');
            
            % Method 2: Ful number of constraints - 2^T constraints
            if T > 16
                demand_M2 = demand_M1;
                generation_M2 = generation_M1;
                SelectionTime_M2 = 0;
                OptimizationTime_M2 = 10^(6);
            else
                [demand_M2,generation_M2,SelectionTime_M2,OptimizationTime_M2] = M2_PowerSetConstraints(T,N,A0(:,:,s),P0(:,s),E0(:,s),Di0(s,:),Gmax,Gmin);
            end
            Obj((n-1)*length(TVec)+t).demand_M2(:,s) = demand_M2;
            Obj((n-1)*length(TVec)+t).generation_M2(:,s) = generation_M2;
            Obj((n-1)*length(TVec)+t).SelectionTime_M2(:,s) = SelectionTime_M2;
            Obj((n-1)*length(TVec)+t).OptimizationTime_M2(:,s) = OptimizationTime_M2;
            fprintf('\nM2 finished! \n');

            % Method 3: Reduced number of constraints - combinatorial search
            [demand_M3,generation_M3,SelectionTime_M3,OptimizationTime_M3] = M3_ReducedConstraints_CombinatorialSearch(T,N,A0(:,:,s),P0(:,s),E0(:,s),Di0(s,:),Gmax,Gmin);
            Obj((n-1)*length(TVec)+t).demand_M3(:,s) = demand_M3;
            Obj((n-1)*length(TVec)+t).generation_M3(:,s) = generation_M3;
            Obj((n-1)*length(TVec)+t).SelectionTime_M3(:,s) = SelectionTime_M3;
            Obj((n-1)*length(TVec)+t).OptimizationTime_M3(:,s) = OptimizationTime_M3;
            fprintf('\nM3 finished! \n');
            
            % Method 4: Reduced number of constraints - greedy search
            [demand_M4,generation_M4,SelectionTime_M4,OptimizationTime_M4] = M4_ReducedConstraints_GreedySearch(T,N,A0(:,:,s),P0(:,s),E0(:,s),Di0(s,:),Gmax,Gmin);
            Obj((n-1)*length(TVec)+t).demand_M4(:,s) = demand_M4;
            Obj((n-1)*length(TVec)+t).generation_M4(:,s) = generation_M4;
            Obj((n-1)*length(TVec)+t).SelectionTime_M4(:,s) = SelectionTime_M4;
            Obj((n-1)*length(TVec)+t).OptimizationTime_M4(:,s) = OptimizationTime_M4;
            fprintf('\nM4 finished! \n');
            
            % Method 5: Reduced number of constraints - greedy search
            [demand_M5,generation_M5,SelectionTime_M5,OptimizationTime_M5] = M5_ReducedConstraints_OneStepSearch(T,N,A0(:,:,s),P0(:,s),E0(:,s),Di0(s,:),Gmax,Gmin);
            Obj((n-1)*length(TVec)+t).demand_M5(:,s) = demand_M5;
            Obj((n-1)*length(TVec)+t).generation_M5(:,s) = generation_M5;
            Obj((n-1)*length(TVec)+t).SelectionTime_M5(:,s) = SelectionTime_M5;
            Obj((n-1)*length(TVec)+t).OptimizationTime_M5(:,s) = OptimizationTime_M5;
            fprintf('\nM5 finished! \n');
        end
        clc;
    end
end

%% Failure rate

for j = 1:length(TVec)*length(NVec)
    for s = 1:NS
        if abs(max(Obj(j).demand_M5(:,s)-Obj(j).demand_M3(:,s))) > 10^(-2) || abs(max(Obj(j).generation_M5(:,s) - Obj(j).generation_M3(:,s))) > 10^(-2)
            Obj(j).FailureNumber_M5 = Obj(j).FailureNumber_M5 + 1;
            Obj(j).Success_Rate_M5 = 1 - Obj(j).FailureNumber_M5/NS;
        end
        
        if abs(max(Obj(j).demand_M4(:,s)-Obj(j).demand_M3(:,s))) > 10^(-2) || abs(max(Obj(j).generation_M4(:,s) - Obj(j).generation_M3(:,s))) > 10^(-2)
            Obj(j).FailureNumber_M4 = Obj(j).FailureNumber_M4 + 1;
            Obj(j).Success_Rate_M4 = 1 - Obj(j).FailureNumber_M4/NS;
        end
    end
end

%% Data to excel -- fixed T
t1 = 1;
SelctionTime_t1 = zeros(6*length(NVec),1);
OptimizationTime_t1 = zeros(6*length(NVec),1);

for n = 1:length(NVec)
    SelctionTime_t1((n-1)*6+1:(n-1)*6+5) = [mean(Obj((n-1)*length(TVec)+t1).SelectionTime_M1); mean(Obj((n-1)*length(TVec)+t1).SelectionTime_M2); mean(Obj((n-1)*length(TVec)+t1).SelectionTime_M3); mean(Obj((n-1)*length(TVec)+t1).SelectionTime_M4); mean(Obj((n-1)*length(TVec)+t1).SelectionTime_M5)];
    OptimizationTime_t1((n-1)*6+1:(n-1)*6+5) = [mean(Obj((n-1)*length(TVec)+t1).OptimizationTime_M1); mean(Obj((n-1)*length(TVec)+t1).OptimizationTime_M2); mean(Obj((n-1)*length(TVec)+t1).OptimizationTime_M3); mean(Obj((n-1)*length(TVec)+t1).OptimizationTime_M4); mean(Obj((n-1)*length(TVec)+t1).OptimizationTime_M5)];
end

t2 = 3;
SelctionTime_t2 = zeros(6*length(NVec),1);
OptimizationTime_t2 = zeros(6*length(NVec),1);

for n = 1:length(NVec)
    SelctionTime_t2((n-1)*6+1:(n-1)*6+5) = [mean(Obj((n-1)*length(TVec)+t2).SelectionTime_M1); mean(Obj((n-1)*length(TVec)+t2).SelectionTime_M2); mean(Obj((n-1)*length(TVec)+t2).SelectionTime_M3); mean(Obj((n-1)*length(TVec)+t2).SelectionTime_M4); mean(Obj((n-1)*length(TVec)+t2).SelectionTime_M5)];
    OptimizationTime_t2((n-1)*6+1:(n-1)*6+5) = [mean(Obj((n-1)*length(TVec)+t2).OptimizationTime_M1); mean(Obj((n-1)*length(TVec)+t2).OptimizationTime_M2); mean(Obj((n-1)*length(TVec)+t2).OptimizationTime_M3); mean(Obj((n-1)*length(TVec)+t2).OptimizationTime_M4); mean(Obj((n-1)*length(TVec)+t2).OptimizationTime_M5)];
end

ExcelDataInput = table(SelctionTime_t1,OptimizationTime_t1,SelctionTime_t2,OptimizationTime_t2);
writetable(ExcelDataInput,'Figures.xlsx','Sheet','FixedTime')

%% Data to excel -- fixed N
n1 = 2;
SelctionTime_n1 = zeros(6*length(TVec),1);
OptimizationTime_n1 = zeros(6*length(TVec),1);

for t = 1:length(TVec)
    SelctionTime_n1((t-1)*6+1:(t-1)*6+5) = [mean(Obj((n1-1)*length(TVec)+t).SelectionTime_M1); mean(Obj((n1-1)*length(TVec)+t).SelectionTime_M2); mean(Obj((n1-1)*length(TVec)+t).SelectionTime_M3); mean(Obj((n1-1)*length(TVec)+t).SelectionTime_M4); mean(Obj((n1-1)*length(TVec)+t).SelectionTime_M5)];
    OptimizationTime_n1((t-1)*6+1:(t-1)*6+5) = [mean(Obj((n1-1)*length(TVec)+t).OptimizationTime_M1); mean(Obj((n1-1)*length(TVec)+t).OptimizationTime_M2); mean(Obj((n1-1)*length(TVec)+t).OptimizationTime_M3); mean(Obj((n1-1)*length(TVec)+t).OptimizationTime_M4); mean(Obj((n1-1)*length(TVec)+t).OptimizationTime_M5)];
end

n2 = 4;
SelctionTime_n2 = zeros(6*length(TVec),1);
OptimizationTime_n2 = zeros(6*length(TVec),1);

for t = 1:length(TVec)
    SelctionTime_n2((t-1)*6+1:(t-1)*6+5) = [mean(Obj((n2-1)*length(TVec)+t).SelectionTime_M1); mean(Obj((n2-1)*length(TVec)+t).SelectionTime_M2); mean(Obj((n2-1)*length(TVec)+t).SelectionTime_M3); mean(Obj((n2-1)*length(TVec)+t).SelectionTime_M4); mean(Obj((n2-1)*length(TVec)+t).SelectionTime_M5)];
    OptimizationTime_n2((t-1)*6+1:(t-1)*6+5) = [mean(Obj((n2-1)*length(TVec)+t).OptimizationTime_M1); mean(Obj((n2-1)*length(TVec)+t).OptimizationTime_M2); mean(Obj((n2-1)*length(TVec)+t).OptimizationTime_M3); mean(Obj((n2-1)*length(TVec)+t).OptimizationTime_M4); mean(Obj((n2-1)*length(TVec)+t).OptimizationTime_M5)];
end

ExcelDataInput = table(SelctionTime_n1,OptimizationTime_n1,SelctionTime_n2,OptimizationTime_n2);
writetable(ExcelDataInput,'Figures.xlsx','Sheet','FixedNumber')

