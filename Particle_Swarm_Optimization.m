clc; clear; close all;

%% DEFINIÇÃO DO PROBLEMA

problem.Function = @(x) pid_separador(x);                       % Valor da Função

problem.nVar = 6;                                               % Numero de Variáveis controladas

problem.VarMax = [10 10 10 10 10 10 ];                          % Limite Superior das Variáveis
problem.VarMin = [-10 -10 -10 -10 -10 -10];                     % Limite Inferior das Variáveis

problem.inic = [-1.7658 -5.0552 0.2224 -0.0999 -0.4839 2.5630]; % Chute Inicial

%% PARÂMETROS DO PSO

params.MaxIt = 30;          % Numero Máximo de Iterações
params.nPop = 10;           % Tamanho da População

params.w = 1;               % Coeficiente Inercial

params.c1 = 2.3;            % Coeficiente de Aceleração Pessoal
params.c2 = 1.8;            % Coeficiente de Aceleração Social

%% CHAMANDO A FUNÇÃO PSO

out = PSO(problem, params);

%% RESULTADOS DA OTIMIZAÇÃO

BestSol = out.BestSol;
Populacao = out.pop;

%% FUNÇÃO PSO

function out = PSO(problem, params)

    %% DEFINIÇÃO DO PROBLEMA

    Function = problem.Function;   % Função Separador

    nVar = problem.nVar;           % Numero de Variáveis

    VarSize = [1 nVar];            % Matriz Numero de Variáveis

    VarMin = problem.VarMin;       % Limite Superior das Variáveis
    VarMax = problem.VarMax;       % Limite Inferior das Variáveis

    %% PARAMETROS DO PSO
    
    MaxIt = params.MaxIt;   % Numero Máximo de Iterações

    nPop = params.nPop;     % Tamanho da População

    w = params.w;           % Coeficiente Inercial 
    
    c1 = params.c1;         % Coeficiente de Aceleração Pessoal
    c2 = params.c2;         % Coeficiente de Aceleração Social

    MaxVelocity = VarMax;
    MinVelocity = VarMin;

    %% INICIALIZANDO VARIÁVEIS

    % Variáveis da Particula
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.BestPosition = [];
    empty_particle.Fitness = [];
    
    % Criando a População
    particle = repmat(empty_particle, nPop, MaxIt);

    % Inicializando Melhor Global e Melhor Fitness
   
    GlobalBest.Position = inf;
    GlobalBest.Fitness = inf;
    

    % Gerando Uma Solução Para Cada Particula
    for i = 1:nPop
        if i == 1
            particle(i).Position = problem.inic;
        else
            particle(i).Position = VarMin(:)' + rand(1,nVar).*(VarMax(:)'-VarMin(:)');
        end
    end
    
    % Inicializando Velocidade
    for i=1:nPop
        particle(i).Velocity = zeros(VarSize);
    end
    
    % Inicializando Melhor Individual
    for i=1:nPop
        
        % Inicializando Fitness
        particle(i).Fitness = Function(particle(i).Position);
        
        % Inicializando Melhor Fitness da Particula   
        if particle(i).Fitness < GlobalBest.Fitness
            GlobalBest.Fitness = particle(i).Fitness;
            GlobalBest.Position = particle(i).Position;
        end
        
        % Inicializando Melhor Posição da Particula
        particle(i).BestPosition = particle(i).Position;        
    end
    
    
    %% INICIANDO PSO

    for it=2:MaxIt
        for i=1:nPop
            
            % Atualizando Velocidade
            particle(i,it).Velocity = w*particle(i,it-1).Velocity ...
                + c1*rand(VarSize).*(particle(i,it-1).BestPosition - particle(i,it-1).Position) ...
                + c2*rand(VarSize).*(GlobalBest.Position - particle(i,it-1).Position);

            % Aplicando os Limites da Velocidade
            for j=1:nVar
                particle(i,it).Velocity(j) = max(particle(i,it).Velocity(j), MinVelocity(j));
                particle(i,it).Velocity(j) = min(particle(i,it).Velocity(j), MaxVelocity(j));
            end

            % Atualizando Posição
            particle(i,it).Position = particle(i,it-1).Position + particle(i,it).Velocity;

            % Aplicando Limites da Posição
            for k=1:nVar
                particle(i,it).Position(k) = max(particle(i,it).Position(k), VarMin(k));
                particle(i,it).Position(k) = min(particle(i,it).Position(k), VarMax(k));
            end
            
            % Evolução
            particle(i,it).Fitness = Function(particle(i,it).Position);

            % Atualizando Melhor Global   
            if particle(i,it).Fitness < GlobalBest.Fitness
                GlobalBest.Fitness = particle(i,it).Fitness;
                GlobalBest.Position = particle(i,it).Position;
            end
            
            % Atualizando Melhor Individual
            if particle(i,it).Fitness < particle(i,it-1).Fitness
                particle(i,it).BestPosition = particle(i,it).Position;
            else
                particle(i,it).BestPosition = particle(i,it-1).BestPosition;
            end

        end
    end
    
    out.pop = particle;
    out.BestSol = GlobalBest;
    
end
