    function maior = Separator(x)

%    x = [-0.2313 -0.1058 3.4650 -3.8772 -7.9190 -6.3831];  % Fitness 1 47.3204    0.0018    0.2906   44.2260
%    x = [-0.0607 -0.0355 0.3292 -1.2523 -3.8693 -1.7440];  % Fitness 1 45.8176    0.0034    1.0640   45.9632
%    x = [-9.9850 -9.9883 9.9212 -9.7620 -9.9698 0.5639];   % Fitness 3 47.6406    0.0014    0.0030   43.4185

    disp(x);

    P1 = x(1);
    I1 = x(2);
    D1 = x(3);
    P2 = x(4);
    I2 = x(5);
    D2 = x(6);


    %% PAR�METROS

    cvg = 120;      % Coeficiente da V�lvula de G�s
    cvl = 1025;     % Coeficiente da V�lvula de L�quido

    p1 = 6;         % Press�o a Jusante da V�lvula de L�quido [bar]
    p2 = 6;         % Press�o a Jusante da V�lvula de G�s [bar]

    rhol = 850;     % Densidade do L�quido [kg/m�]
    rhoh = 999.19;  % Densidade da �gua [kg/m�]

    g = 9.81;       % Gravidade [m/s�]

    v = 56.5487;    % Volume do Separador [m�]
    d = 3;          % Di�metro do Separador [m]
    c = 8;          % Comprimento do Separador [m]

    MMar = 0.02897; % Massa Molar do Ar [kg/mol]
    MMg = 0.021;    % Massa Molar do G�s [kg/mol]

    %% CONDI��ES ESTADO ESTACION�RIO

    lin = 0.165;      % Vaz�o de L�quido na Entrada do Separador [m�/s]
    lout = 0.165;     % Vaz�o de L�quido na Sa�da do Separador [m�/s]

    hl = 2;           % Altura de L�quido no Separador [m]
    p = 8;            % Press�o de G�s no Separador [bar]

    gin = 0.1;        % Vaz�o de G�s na Entrada do Separador [m�/s]
    gout = 0.1;       % Vaz�o de G�s na Sa�da do Separador [m�/s]

    xl = 0.5;         % Abertura da V�lvula de L�quido
    xg = 0.5;         % Abertura da V�lvula de G�s

    vl = 40.0483;     % Volume de L�quido no Separador [m�]

    t = 303.15;       % Temperatura [K]
    R = 8.3144621e-5; % Constante Universal dos G�ses [ bar m�/(gmol.K) ]


    %% Fun��es de Transfer�ncia

    % ----------------------------- ALTURA ------------------------------------

    % Declarando as Vari�veis
    syms LIN LOUT HL XL P

    % Equa��es de Conserva��o
    Ff  = (LIN-LOUT)/(2*c*sqrt(HL*(d-HL)));                                 % Balan�o de Massa para a Fase L�quida - dh/dt = f
    Flo = 2.4e-4*XL*cvl*sqrt((P+(rhol*g*HL*10^-5)-p1)/(rhol/rhoh));         % Equa��o da V�lvula de L�quido - Lout

    % DERIVADAS - Altura (serao utilizadas para calcular as FTs)

    Ff_li = subs(diff(Ff,LIN),{LIN,LOUT,HL},[lin,lout,hl]);                 % df/dLin
    Ff_lo = subs(diff(Ff,LOUT),{LIN,LOUT,HL},[lin,lout,hl]);                % df/dLout
    Ff_hl = subs(diff(Ff,HL),{LIN,LOUT,HL},[lin,lout,hl]);                  % df/dhl

    Flo_xl = sym2poly(subs(diff(Flo,XL),{XL,P,HL},[xl,p,hl]));              % dLout/dxl
    Flo_p  = sym2poly(subs(diff(Flo,P),{XL,P,HL},[xl,p,hl]));               % dLout/dp
    Flo_hl = subs(diff(Flo,HL),{XL,P,HL},[xl,p,hl]);                        % dLout/dhl

    % Criando as Fun��es de Transfer�ncia  a = 1/tau
    a = -sym2poly((Ff_hl+(Ff_lo*Flo_hl)));

    % FT 1 - ALTURA DO L�QUIDO x FRA��O DE ABERTURA DA V�LVULA DE L�QUIDO  Hl(s)/Xl(s)
    k1   = sym2poly(Ff_lo * Flo_xl)/a;
    tau1 = 1/a;
    g1   = tf(k1,[tau1 1]);

    % FT 2 - ALTURA DO L�QUIDO x VAZ�O DE ENTRADA DE L�QUIDO               Hl(s)/Lin(s)
    k2   = sym2poly(Ff_li)/a;
    tau2 = 1/a;
    g2   = tf([k2],[tau2 1]);

    % FT 3 ALTURA DE L�QUIDO x PRESS�O                                     Hl(s)/P(s)
    k3   = sym2poly(Ff_lo * Flo_p)/a;
    tau3 = tau1;
    g3   = tf(k3,[tau3 1]);


    % ----------------------------- PRESS�O -----------------------------------

    % Declarando Vari�veis
    syms GIN GOUT VL XG T

    % Equa��es de Conserva��o
    Fvl   = (c*(d^2)/4)*(acos((d-2*HL)/d)-((2*sqrt((d-HL)*HL)/d)*(((d-(2*HL))/d))));    % Volume de L�quido em fun��o de h(t) - Vl
    Fg    = P*(GIN-GOUT+LIN-LOUT)/(v-VL);                                               % Balan�o de Massa da Fase Gasosa - dp/dt = g
    Fgout = 2.4e-4*XG*cvg*sqrt((P-p2)*(P+p2)*T*MMar/MMg/P^2);                           % Equa��o da V�lvula de G�s - Gout

    % DERIVADAS - Press�o (serao utilizadas para calcular as FTs)

    Fvl_hl = subs(diff(Fvl,HL),HL,hl);                                                  % dVl/dt

    Fg_li = subs(diff(Fg,LIN),{GIN,GOUT,LIN,LOUT,VL,P},[gin,gout,lin,lout,vl,p]);   % dg/dLin
    Fg_p  = subs(diff(Fg,P),{GIN,GOUT,LIN,LOUT,VL,P},[gin,gout,lin,lout,vl,p]);     % dg/dp
    Fg_lo = subs(diff(Fg,LOUT),{GIN,GOUT,LIN,LOUT,VL,P},[gin,gout,lin,lout,vl,p]);  % dg/dLout
    Fg_gi = subs(diff(Fg,GIN),{GIN,GOUT,LIN,LOUT,VL,P},[gin,gout,lin,lout,vl,p]);   % dg/dGin
    Fg_go = subs(diff(Fg,GOUT),{GIN,GOUT,LIN,LOUT,VL,P},[gin,gout,lin,lout,vl,p]);  % dg/dGout
    Fg_vl = subs(diff(Fg,VL),{GIN,GOUT,LIN,LOUT,VL,P},[gin,gout,lin,lout,vl,p]);    % dg/dVl

    Fgo_p  = subs(diff(Fgout,P),{XG,P,T},[xg,p,t]);                                     % dGout/dp
    Fgo_xg = subs(diff(Fgout,XG),{XG,P,T},[xg,p,t]);                                    % dGout/dxg
    Fgo_t  = subs(diff(Fgout,T),{XG,P,T},[xg,p,t]);                                     % dGout/dT

    % Criando as Fun��es de Transfer�ncia  den = 1/tau
    den = -sym2poly( Fg_p + (Fg_lo*Flo_p) + (Fg_go*Fgo_p) );

    % FT 4 - PRESS�O x VAZ�O DE ENTRADA DE L�QUIDO                        P(s)/Lin(s)
    k4   = sym2poly(Fg_li) /den;
    tau4 = 1 /den;
    g4   = tf(k4,[tau4 1]);

    % FT 5 - PRESS�O x VAZ�O DE ENTRADA DE G�S                             P(s)/Gin(s)
    k5   = sym2poly(Fg_gi)/den;
    tau5 = tau4;
    g5   = tf(k5,[tau5 1]);

    % FT 6 - PRESS�O x FRA��O DE ABERTURA DA V�LVULA DE L�QUIDO            P(s)/Xl(s)
    k6   = sym2poly(Fg_lo* Flo_xl) /den;
    tau6 = tau4;
    g6   = tf(k6,[tau6 1]);

    % FT 7 - PRESS�O x ALTURA DO L�QUIDO                                   P(s)/Hl(s)
    k7   = sym2poly((Fg_lo* Flo_hl)+(Fg_vl*Fvl_hl))/den;
    tau7 = tau4;
    g7   = tf(k7,[tau7 1]);

    % FT 8 - PRESS�O x FRA��O DE ABERTURA DE G�S                           P(s)/Xg(s)
    k8   = sym2poly(Fg_go* Fgo_xg) /den;
    tau8 = tau4;
    g8   = tf(k8,[tau8 1]);

    % FT 9 - PRESS�O x TEMPERATURA                                         P(s)/T(s)
    k9   = sym2poly(Fg_go* Fgo_t) /den;
    tau9 = tau4;
    g9   = tf(k9,[tau9 1]);


    %% SIMULINK

    % Definindo a Fun��o do SIMULINK
    simNome = 'Separador';

    % Tempo simula��o
    simTime = 30000;

    % para ler as variaveis da funcao e n so do workspace
    old_opt = simget(simNome);
    opt = simset(old_opt,'SrcWorkspace','current');

    % Rodar o Modelo do SIMULINK
    try
        sim(simNome,[],opt);
        disp('try');
    catch
        maior = Inf;
        disp('catch');
        return
    end

    % Dados de Sa�da do SIMULINK
    tt      = P_smlk.time;            % Tanto faz usar P ou H, o tempo � o mesmo
    h_final = H_smlk.data;            % Varia��o do N�vel
    p_final = P_smlk.data;            % Varia��o da Press�o
    x_l     = Xl_smlk.data;           % Varia��o da Abertura da V�lvula de L�quido
    x_g     = Xg_smlk.data;           % Varia��o da Abertura da V�lvula de G�s


    mOut = [h_final p_final];
    xOut = [x_l x_g];

%% AN�LISE DAS SA�DAS E DETERMINA��O DO FITNESS

    compr = length(tt);                             % Comprimento Vetor Reposta
    faixa = 7*round(compr/8);                       % Analisar Parte Est�vel (�ltimo 1/8)
    
    aux = ((p_final+p)-p1)/(rhol/rhoh);             % Vari�vel para o c�lculo da vaz�o de sa�da
    loutt = 2.4*10^(-4)*(xl+x_l)*cvl.*sqrt(aux);    % Vaz�o de Sa�da de L�quido
    
    % Sa�da de L�quido
    max_lout_err = 100*max(abs(loutt-lout)/lout);   % Erro Percentual M�ximo do Vetor de Desvio de Sa�da de L�quido
    
    % Abertura das V�lvulas 
    analisev    = xOut(faixa:end,:);                % Valores do Desvio das V�lvulas (na regi�o est�vel)
    maxv        = max(abs(analisev));               % M�ximo Desvio Absoluto de Cada V�lvula 
    max_val_err = 100*max(maxv)/xl;                 % Erro Percentual M�ximo do Vetor de Desvio de Abertura das V�lvulas
    
    % Press�o e N�vel
    analises  = mOut(faixa:end,:);                  % Valores do Desvio do N�vel e da Press�o (na regi�o est�vel)
    maxs      = max(abs(analises));                 % Maximo Desvio Absoluto de Cada Sa�dal
    max_h_err = 100*max(maxs(:,1)/hl);              % Erro Percentual M�ximo do Vetor de Desvio de N�vel
    max_p_err = 100*max(maxs(:,2)/p);               % Erro Percentual M�ximo do Vetor de Desvio de Press�o

    % FITNESS
%    maior = max([max_val_err max_p_err max_h_err max_lout_err]);                % Maior Valor entre os Erros Percentuais M�ximos
%    maior = max_val_err + max_p_err + max_h_err + max_lout_err;                 % Soma dos Erros
    maior = sqrt(max_val_err^2 + max_p_err^2 + max_h_err^2 + max_lout_err^2);   % Raiz da Soma dos Quadrados dos Erros

    
    disp(maior);              
end