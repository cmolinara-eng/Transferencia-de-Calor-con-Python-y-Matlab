    
    %% Ejemplo 4.1.
    clear;
    % Parámetros
    A_s = 20;                % [m^2]
    esp = 0.30;              % [m]
    T_in = 25 + 273.15;      % T[K]
    k = 1;                   % [W/m·K]
    
    % Rango de temperaturas exteriores (invierno a verano)
    T_inv = -15 + 273.15;    %  [K]
    T_ver = 38 + 273.15;     %  [K]
    
    T_ev = linspace(T_inv, T_ver, 100);
    
    Q_a = - A_s * k * (T_in - T_ev) / esp;  % [W]
    
    % Gráfico
    figure('Color', 'w');
    hold on;
    plot(T_ev - 273.15, Q_a, 'k--', 'LineWidth', 2, 'DisplayName', sprintf('Base (k=%.2f)', k));
    
    % Conductividades adicionales
    k_materiales = [0.75, 1.25, 3.0, 5];
    
    
    for i = 1:length(k_materiales)
        k_val = k_materiales(i);
        Q_b = - A_s * k_val * (T_in - T_ev) / esp;
        plot(T_ev - 273.15, Q_b, 'LineWidth', 1.5, 'DisplayName', sprintf('(k=%.2f)', k_val));
    end
    
    xlabel('Temperatura exterior (°C)', 'FontSize', 12);
    ylabel('Flujo de calor (W)', 'FontSize', 12);
    title('Flujo de calor en función de la temperatura exterior', 'FontSize', 14);
    grid off;
    yline(0, '--', 'Color', [0.5 0.5 0.5]);
    legend('Location', 'best');
    hold off;



    %% Ejemplo 4.2.
    % Datos
    clear;
    T_inf = 45 + 273.15;  % [K]
    T_s = 280 + 273.15;   % [K]
    V = [1 2 4 8 12];     % [m/s]
    D = 25e-3;            % [m]
    Q = [310 470 720 1150 1600]; % [W/m]
    
    % Cálculo del coeficiente convectivo
    h = Q ./ (pi * D * (T_s - T_inf));
    
    % Ajuste del modelo h = C * V^n
    modelo = @(p, V) p(1) * V.^p(2);
    p0 = [1 1];
    p = lsqcurvefit(modelo, p0, V, h);
    C = p(1);
    n = p(2);
    
    % R²
    h_pred = modelo(p, V);
    res = h - h_pred;
    R2 = 1 - sum(res.^2) / sum((h - mean(h)).^2);
    
    % Resultados
    fprintf('C = %.4f\n', C);
    fprintf('n = %.4f\n', n);
    fprintf('R² = %.4f\n', R2);
    
    % Gráfico
    V_fit = linspace(min(V), max(V), 100);
    h_fit = modelo(p, V_fit);
    
    figure('Color', 'w');
    scatter(V, h, 'filled'); hold on;
    plot(V_fit, h_fit, 'r', 'LineWidth', 1.5);
    xlabel('Velocidad (V) [m/s]');
    ylabel('Coeficiente h [W/m^2·K]');
    title('h = C·V^n');
    legend('Datos', sprintf('Ajuste: h = %.2f·V^{%.2f}', C, n), 'Location', 'best');
    grid off;
    hold off;

    %% Ejemplo 4.3.
    clear;
    % Planteamos los parametros
    D = 0.15;                     % Diametro [m]
    A = pi * D^2;            % Surface area of sphere [m^2]
    sigma = 5.67e-8;             % Constante de Stefan-Boltzmann [W/m^2.K^4]
    T_alr = 10+273;                  % Temperatura alrededor [K]
    
    % Establecemos el limite de temperatura para los paquetes
    T_C = 35:1:95;               % [°C]
    T_K = T_C + 273.15;          % Conversión [K]
    
    % Emisividades 
    emisividades = 0.15:0.02:0.35;
    
    % Inicializamos la grafica
    figure;
    hold on;
    
    % Bucle para iterar sobre los valores de emisividad
    for eps = emisividades
        Q = eps * sigma * A * (T_K.^4 - T_alr^4); % [W]
        plot(T_C, Q, 'DisplayName', sprintf('\\epsilon = %.2f', eps));
    end
    
    % Plot 
    xlabel('Temperatura de la superficie (°C)');
    ylabel('Disipación (W)');
    title('Disipación de calor para diferentes ε');
    legend('Location', 'northwest');
    grid off;




    %% Ejemplo 4.4.
    clear; clc; close all;
    % Datos
    h = 25;                    % W/m^2·K
    T_inf = 25 + 273.15;       % K
    T_alr = 25 + 273.15;       % K
    T_s_int = 800;             % K
    eps = 0.8;                 % Emisividad
    r1 = 0.12 / 2;             % m
    boltz = 5.67e-8;           % W/m^2·K^4
    esp_max = 0.14;            % m
    k = 0.09;                  % W/m·K
    
    % Función a resolver
    minfunc = @(T_s2, esp) ...
        ( h * (T_s2 - T_inf) + ...                          % Convección
          eps * boltz * (T_s2.^4 - T_alr.^4) - ...          % Radiación
          k * (T_s_int - T_s2) ./ ((r1 + esp) .* log((r1 + esp) ./ r1)) ); % Conducción
    
    
    % Vector de espesores
    espesores = linspace(0.001, esp_max, 50);  % Evitar log(0)
    T2_valores = zeros(size(espesores));

    % Solución con fsolve
    options = optimset('Display','off');  % silenciar salida
    
    for j = 1:length(espesores)
        esp = espesores(j);
        if j == 1
            v_in = T_s_int - 400;
        else
            v_in = T2_valores(j-1);
        end
        % Resolver ecuación
        T2_valores(j) = fsolve(@(T) minfunc(T, esp), v_in, options);
    end
    
    % Cálculo de flujo de calor por unidad de longitud
    r2 = r1 + espesores;
    q_cond = (k * 2 * pi * (T_s_int - T2_valores)) ./ log(r2 ./ r1); % W/m
    
    % Coordenada adimensional
    cord_ad = espesores ./ esp_max;

    figure;
    plot(cord_ad, q_cond, 'o-', 'Color', [0 0.5 0.5], 'MarkerFaceColor', [0 0.5 0.5]);
    title('Calor vs. Espesor del aislamiento');
    xlabel('Coordenada adimensional');
    ylabel('Calor (W/m)');
    grid off;
    set(gca, 'FontSize', 11);
    
    figure;
    plot(cord_ad, T2_valores, 'o-', 'Color', [0 0.5 0.5], 'MarkerFaceColor', [0 0.5 0.5]);
    title('Temperatura de la superficie externa vs. Espesor del aislamiento');
    xlabel('Coordenada adimensional');
    ylabel('Temperatura superficie externa T_{s2} [K]');
    grid off;
    set(gca, 'FontSize', 11);


    %% Ejemplo 4.5.
    clear;
    %Placa plana
    K = 16.3;           % Conductividad [W/(m*K)]
    T1 = 473;           % [K]
    T2 = 373;           % [K]
    L = 0.80;           % [m]
    x1 = 0;
    x2 = L;
    % Radios
    r1 = 0.10;          % [m]
    r2 = 0.025;         % [m]
    % Coeficientes del radio lineal r(x) = mx + b
    m = (r2 - r1)/L;
    b = r1;
    % Calor transmitido
    n = -m * pi * K * (T2 - T1);
    d = (1 / (m * x1 + b)) - (1 / (m * x2 + b));
    Qx = n / d;   % Resultado en [W]
    % Temperatura en x = 0.5 m
    x_eval = 0.50;
    term1 = (1 / (m * L + b)) - (1 / b);
    term2 = (1 / (m * x_eval + b)) - (1 / b);
    T = ((T2 - T1) / term1) * term2 + T1;  % [K]
    % Perfil de temperatura
    x_vals = linspace(0, L, 100);   % Rango de x de 0 a 0.8 m
    T_vals = ((T2 - T1) ./ ((1./(m*L + b)) - (1./b))) .* ((1./(m*x_vals + b)) - (1./b)) + T1;
    fprintf('Calor transmitido Qx = %.1f W\n', Qx);
    fprintf('Temperatura en x = 0.5 m: T = %.2f K\n', T);
    % Perfil de temperatura
    figure;
    plot(x_vals, T_vals, 'r', 'LineWidth', 2);
    xlabel('Longitud x [m]');
    ylabel('Temperatura [K]');
    title('Perfil de temperatura');
    grid off;

    %% Ejemplo 4.6.
    clear;
    %Cilindro
    % Parámetros
    R = 0.004;              % [m]
    k = 16;                 % [W/m K]
    Phi = 5e6;              % Generación de calor [W/m^3]
    T_s = 280;              % [K]
    h = 1;                  % Altura del cilindro [m]
    % Vector de radios
    r = linspace(0, R, 200);
    % Perfil de temperatura
    T_r = T_s + (Phi * R^2 / (4 * k)) * (1 - (r / R).^2);
    % Flujo de calor
    q_r = pi * h * Phi * r.^2;
    % perfil de temperatura
    figure;
    plot(r, T_r, 'b', 'LineWidth', 2);
    xlabel('r (m)');
    ylabel('T (K)');
    title('Perfil de temperatura');
    grid on;
    % Flujo de calor
    figure;
    plot(r, q_r, 'b', 'LineWidth', 2);
    xlabel('r (m)');
    ylabel('Q_r (W/m²)');
    title('Flujo de calor radial');
    grid off;

    %% Ejemplo 4.7 
    clc; clear; close all;
    
    % Datos
    R = 0.1;               % [m]
    Gen = 5e4;             % [W/m^3]
    k = 25;                % [W/m K]
    T_inf = 350;           % [K]
    h = 100;               % [W/m^2*K]
    
    
    radios = linspace(1e-6, R, 100); 
    
    % Solución Analítica
    T_eval = T_inf + Gen * (R^2 - radios.^2) / (6 * k) + Gen * R / (3 * h);
    
    figure('Name', 'Perfil de temperatura - Analítico', 'NumberTitle', 'off');
    plot(radios, T_eval, 'b-', 'LineWidth', 2);
    xlabel('r [m]');
    ylabel('T [K]');
    title('Perfil de temperatura (Analítica)');
    legend('Analítica');
    grid on;
    set(gca, 'FontSize', 10);
    
    
    r_sol = linspace(1e-6, R, 100); 
    
    
    y_init = [T_inf; 0]; 
    solinit = bvpinit(r_sol, y_init);
    % ------------------------------------------
    
    % Definición de funciones
    odefun = @(r, y) odefun_esfera(r, y, Gen, k);
    bcfun_handle = @(ya, yb) bcfun_esfera(ya, yb, k, h, T_inf);
    
    % Solución
    sol = bvp4c(odefun, bcfun_handle, solinit);
    
    % Evaluación
    T_num_eval = deval(sol, radios);
    T_num = T_num_eval(1,:); 
    
    figure('Name', 'Comparación de Soluciones', 'NumberTitle', 'off');
    plot(radios, T_num, 'r-', 'LineWidth', 2, 'DisplayName', 'Numérica (bvp4c)');
    hold on;
    plot(radios, T_eval, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Analítica');
    hold off;
    xlabel('r [m]');
    ylabel('T [K]');
    title('Comparación: Analítica vs Numérica (bvp4c)');
    legend('show', 'Location', 'best');
    grid on;
    set(gca, 'FontSize', 10);


    %% Ejemplo 4.8.
    clear;
     % Parametros
    Di = 0.05;          % [m]
    e_tub = 0.003;      % [m]
    k_tub = 12.8;        % [W/m K]
    k_aislante = 0.23;   % [W/m K]
    h_i = 433;           % Coef. convección interno [W/m^2 K]
    h_o = 8;             % Coef. convección externo [W/m^2 K]
    T_i = 7;             % [°C]
    T_o = 25;            % [°C]
    e_aislante = 0.012;  % [m]
    r1 = Di / 2;                  % Radio interior del tubo
    r2 = r1 + e_tub;              % Radio exterior del tubo
    r3 = r2 + e_aislante;         % Radio exterior del aislante
    % Sin aislante
    
    R_conv_i = 1 / (2*pi*r1*h_i);
    R_cond_acero = log(r2/r1) / (2*pi*k_tub);
    R_conv_o = 1 / (2*pi*r2*h_o);
    R_total_sin = R_conv_i + R_cond_acero + R_conv_o;
    % Flujo de calor sin aislamiento
    q_sin = (T_o - T_i) / R_total_sin;
    % Con aislante
    R_cond_aislante = log(r3/r2) / (2*pi*k_aislante);
    R_conv_o_aislado = 1 / (2*pi*r3*h_o);
    R_total_con = R_conv_i + R_cond_acero + R_cond_aislante + R_conv_o_aislado;
    q_con = (T_o - T_i) / R_total_con;
    % Radio crítico
    r_critico = k_aislante / h_o;
    fprintf('--- Resultados ---\n');
    fprintf('\n Sin aislamiento:\n');
    fprintf('R_conv_i  = %.4e K/W·m\n', R_conv_i);
    fprintf('R_cond    = %.4e K/W·m\n', R_cond_acero);
    fprintf('R_conv_o  = %.4e K/W·m\n', R_conv_o);
    fprintf('R_total   = %.2f K/W·m\n', R_total_sin);
    fprintf('q_sin     = %.2f W/m\n', q_sin);
    fprintf('\n Con aislamiento de 10 mm:\n');
    fprintf('R_cond_aislante = %.2f K/W·m\n', R_cond_aislante);
    fprintf('R_conv_o_aislado = %.2f K/W·m\n', R_conv_o_aislado);
    fprintf('R_total = %.2f K/W·m\n', R_total_con);
    fprintf('q_con   = %.2f W/m\n', q_con);
    fprintf('\nRadio crítico total (tubo + aislante): %.4f m (%.1f mm)\n', r_critico, r_critico*1000);
    fprintf('\nEspesor crítico del aislante: %.4f m (%.1f mm)\n', (r_critico - r2), (r_critico - r2)*1000);
    % q vs espesor del aislante
    
    espesores_m = linspace(0, 0.03, 300);  % de 0 a 30 mm
    q_vs_e = zeros(size(espesores_m));    % inicializamos
    for i = 1:length(espesores_m)
       r3_var = r2 + espesores_m(i);      % radio exterior variable
       R_cond_aislante = log(r3_var / r2) / (2 * pi * k_aislante);
       R_conv_o_aislado = 1 / (2 * pi * r3_var * h_o);
       R_total_var = R_conv_i + R_cond_acero + R_cond_aislante + R_conv_o_aislado;
       q_vs_e(i) = (T_o - T_i) / R_total_var;
    end
    figure;
    plot(espesores_m*1000, q_vs_e, 'b-', 'LineWidth', 2); hold on;
    xline((r_critico - r2)*1000, '--r', 'Espesor crítico');
    xlabel('Espesor del aislante [mm]');
    ylabel('Flujo de calor q [W/m]');
    title('Espesor del aislamiento vs flujo de calor');
    grid off;

    %% Ejemplo 4.9
    clc; clear; close all;
    
    % Parámetros
    params_fin.P = 0.066;
    params_fin.A_c = 2e-4;
    params_fin.L = 0.1;
    params_fin.k = 150;
    params_fin.w = 0.025;
    params_fin.t = 0.008;
    params_fin.D = 0.0; 
    
    % FluidParams
    params_fluid.h = 35;
    params_fluid.T_inf = 300;

    casos = {"longitud infinita", "tconocida", "adiabática", "general"};
    L = params_fin.L;
    T_b = 500;
    X_vals = linspace(0, L, 100);
    
    figure('Name', 'Distribución de Temperatura en la Aleta', 'NumberTitle', 'off', 'Position', [100 100 800 500]);
    hold on;
    for i = 1:length(casos)
        caso = casos{i};
        T_x_list = zeros(size(X_vals));
        
        % Inicializamos res fuera del bucle para tener el último resultado después de iterar
        res = struct('q_f', 0, 'Temperatura', 0, 'eta_A', nan, 'epsilon', nan, 'Caso_efectivo', '');
    
        for j = 1:length(X_vals)
            X = X_vals(j);
            % Llama a la función local
            res = superficies_q_th_matlab(caso, params_fin, params_fluid, T_b, X);
            T_x_list(j) = res.Temperatura;
        end
        
        % Preparar leyenda con el último resultado de la iteración (res)
        q_f = res.q_f;
        eta_A = res.eta_A;
        epsilon = res.epsilon;
        caso_efectivo = res.Caso_efectivo;
    
        if ~isnan(eta_A)
            % Uso de \eta para el símbolo de eficiencia en la leyenda
            lbl = sprintf('%s | q=%.2f W, \\eta=%.3f, \\epsilon=%.3f', caso_efectivo, q_f, eta_A, epsilon);
        else
            lbl = sprintf('%s | q=%.2f W', caso_efectivo, q_f);
        end
    
        plot(X_vals, T_x_list, 'LineWidth', 2, 'DisplayName', lbl);
    end
    plot([0 L], [params_fluid.T_inf params_fluid.T_inf], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'DisplayName', 'T_{\infty}');
    xlabel('x [m]');
    ylabel('Temperatura [K]');
    title('Distribución de temperatura en la aleta rectangular');
    legend('show', 'FontSize', 9);
    grid on;
    grid minor;
    box on;
    hold off;

    function res = superficies_q_th_matlab(Caso, Params_fin, params_fluid, T_b, X)
    
        C = lower(Caso);
    
        % Parámetros de la Aleta
        P = Params_fin.P;
        A_c = Params_fin.A_c;
        L = Params_fin.L;
        k = Params_fin.k;
        
        % Parámetros del Fluido
        h = params_fluid.h;
        T_inf = params_fluid.T_inf;
    
        % Parámetros Comunes
        theta_b = T_b - T_inf;
        m = sqrt(h * P / (k * A_c));
        M = sqrt(h * P * k * A_c) * theta_b;
    
        % Valores predeterminados
        q_f = 0;
        theta_ratio = 0;
        caso_efectivo = '';
    
        % --- Cuatro Casos ---
        if strcmp(C, 'longitud infinita')
            q_f = M;
            theta_ratio = exp(-m*X);
            caso_efectivo = "longitud infinita";
    
        elseif strcmp(C, 'tconocida')
            q_f = h * P * L * theta_b; 
            theta_ratio = 1 - X/L;
            caso_efectivo = "T punta conocida (Aprox. lineal)";
    
        elseif strcmp(C, 'adiabática')
            q_f = M * tanh(m * L);
            theta_ratio = cosh(m * (L - X)) / cosh(m * L);
            caso_efectivo = "punta adiabática";
    
        elseif strcmp(C, 'general')
            h_mk = h/(m*k);
            num_q = sinh(m*L) + h_mk * cosh(m*L);
            den_q = cosh(m*L) + h_mk * sinh(m*L);
            q_f = M * num_q / den_q;
            
            num_theta = cosh(m*(L-X)) + h_mk * sinh(m*(L-X));
            den_theta = cosh(m*L) + h_mk * sinh(m*L);
            theta_ratio = num_theta / den_theta;
            caso_efectivo = "condición general";
    
        else
            error('Caso inválido: %s', Caso);
        end
    
        T_x = T_inf + theta_b * theta_ratio;
    
        % --- Eficiencia y Efectividad (Solo en Adiabática y General) ---
        eta_A = nan;
        epsilon = nan;
    
        if strcmp(C, "adiabática") || strcmp(C, "general")
            A_s = P*L + Params_fin.A_c; % Uso de A_c del struct Params_fin
            A_a = Params_fin.A_c;
            eta_A = q_f / (h*A_s*theta_b);
            epsilon = q_f / (h*A_a*theta_b);
        end
    
        % Ensamblar resultado como estructura
        res = struct( ...
            'q_f', q_f, ...
            'Temperatura', T_x, ...
            'eta_A', eta_A, ...
            'epsilon', epsilon, ...
            'Caso_efectivo', caso_efectivo ...
        );
    end


%% --- FUNCIONES AUXILIARES 

% Función diferencial para la esfera (Ejemplo 4.7)
function dydr = odefun_esfera(r, y, Gen, k)
    % y(1) es T (Temperatura)
    % y(2) es dT/dr (Pendiente)
    dydr = zeros(2,1);
    dydr(1) = y(2);
    
    if r == 0
        % Tratamiento de la singularidad en el centro (L'Hopital)
        dydr(2) = -Gen / (3*k); 
    else
        % Ecuación esférica estándar: T'' + (2/r)T' + Gen/k = 0
        dydr(2) = -2 * y(2) / r - Gen / k;
    end
end

% Condiciones de frontera (Ejemplo 4.7)
function res = bcfun_esfera(ya, yb, k, h, T_inf)
    % ya -> condiciones en el inicio (r cerca de 0)
    % yb -> condiciones en el final (r = R)
    
    % 1. Simetría en el centro
    cond_centro = ya(2); 
    
    % 2. Convección en la superficie
    cond_superficie = -k * yb(2) - h * (yb(1) - T_inf);
    
    res = [cond_centro; cond_superficie];
end







