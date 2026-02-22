% Ejemplo 7.1

    clear; clc; close all;
    % Datos del problema
    D = 0.20;                 %[m]
    Ts = 500;                 %[K]
    Tinf = 310;               %[K]
    Tf = (Ts + Tinf) / 2;     %[K]
    
    % Propiedades del aire a Tf = 400 K
    k = 0.0331;               % [W/mK]
    nu = 25.5e-6;             % v. cinematica[m²/s]
    mu = 22.52e-6 ;           % v. dinamica [kg/m s] 
    cp = 1009;                % [J/kg K]
    beta = 1 / Tf;            % [1/K]
    g = 9.81;                 % [m/s²]
    
    Pr = (mu * cp)/k;
    
    % Número de Rayleigh
    Ra_D = (g * beta * (Ts - Tinf) * D^3) / (nu^2) * Pr;
    
    % Correlación (ecuación 80)
    term = (0.559 / Pr)^(9/16);
    Nu_D = 0.36 + (0.518 * Ra_D^(1/4)) / ((1 + term)^(4/9));
    
    % Cálculo del coeficiente convectivo promedio
    h_c = (k / D) * Nu_D;     % [W/m²K]
    
    % Pérdida de calor por unidad de longitud de tubo
    A = pi * D;
    Q = h_c * A * (Ts - Tinf);  % [W/m]
    
    % Resultados
    fprintf('Resultados:\n');
    fprintf('Temperatura de película: Tf = %.1f K\n', Tf);
    fprintf('Número de Prandalt: Pr = %.2e\n', Pr);
    fprintf('Número de Rayleigh: Ra_D = %.2e\n', Ra_D);
    fprintf('Número de Nusselt: Nu_D = %.2f\n', Nu_D);
    fprintf('Coeficiente convectivo: h_c = %.2f W/m² K\n', h_c);
    fprintf('Pérdida de calor por unidad de longitud: Qdot = %.2f W/m\n', Q);

    % -----------------------------
    % Parámetros
    % -----------------------------
    Pr = 0.7;                              % Prandtl típico del aire
    Ra = logspace(-5,12,400);              % Vector Rayleigh
    
    % -----------------------------
    % FACTORES AUXILIARES
    % -----------------------------
    Psi_vert = (1 + (0.492/Pr)^(9/16))^(-16/9);
    Psi_cil  = (1 + (0.559/Pr)^(9/16))^(-4/9);
    Psi_esf  = (1 + (0.469/Pr)^(9/16))^(-4/9);
    
    % -----------------------------
    % PARED VERTICAL
    % -----------------------------
    Nu_vert = zeros(size(Ra));
    
    for i = 1:length(Ra)
        if Ra(i) < 1e9
            Nu_vert(i) = 0.68 + 0.670*(Ra(i)*Psi_vert)^(1/4);
        else
            Nu_vert(i) = 0.68 + 0.67*(Ra(i)*Psi_vert)^(1/4) * ...
                         (1 + (1.6e-8*Ra(i)*Psi_vert)^(1/12));
        end
    end
    
    % -----------------------------
    % CILINDRO HORIZONTAL
    % -----------------------------
    Nu_cil = zeros(size(Ra));
    
    for i = 1:length(Ra)
        if Ra(i) <= 1e9
            Nu_cil(i) = 0.36 + (0.518*Ra(i)^(1/4))/((1+(0.559/Pr)^(9/16))^(4/9));
        else
            Nu_cil(i) = (0.60 + 0.387*(Ra(i)/((1+(0.559/Pr)^(9/16))^(16/9)))^(1/6))^2;
        end
    end
    
    % -----------------------------
    % ESFERA
    % -----------------------------
    Nu_esf = 2 + (0.589*Ra.^(1/4))./((1+(0.469/Pr)^(9/16)).^(4/9));
    
    % -----------------------------
    % GRÁFICA
    % -----------------------------
    figure
    loglog(Ra,Nu_vert,'b','LineWidth',2); hold on
    loglog(Ra,Nu_cil,'r','LineWidth',2)
    loglog(Ra,Nu_esf,'Color',[1 0.6 0],'LineWidth',2)
    
    xlabel('Ra','FontSize',12)
    ylabel('Nu','FontSize',12)
    legend('Pared vertical','Cilindro horizontal','Esfera','Location','southeast')
    
    set(gca,'FontSize',11)