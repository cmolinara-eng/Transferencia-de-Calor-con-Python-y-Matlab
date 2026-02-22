    %% Ejemplo 5.1.
    % Pelicula de mercurio
    clear; clc;
    % Datos
    Ts  = 75;             % [C]
    T0  = 25;             % [C]
    e   = 2e-3;           % [m]
    b   = 5;              % [m]
    L   = 3;              % [m]
    rho = 13412;          % [kg/m^3]
    k   = 9.16;           % [W/m-K]
    Cp  = 137.8;          % [J/kg-K]
    Qv  = 5/3600;         % [m^3/s]  (5 m^3/h)
    
    % Propiedades efectivas
    alpha   = k/(rho*Cp);           % [m^2/s]
    Asec    = e*b;                  % [m^2]
    v_med   = Qv/Asec;              % [m/s]
    v_max   = 1.5*v_med;            % perfil laminar
    
    % Funcion temperatura cerrada
    Tfun = @(X,Z) T0 + (1 - erf( X./sqrt(4*alpha*Z./v_max) )) .* (Ts - T0);
    
    % Temperaturas pedidas
    Zs = [1,2,3];             % [m]
    Xs = [0, e/2, e];         % pared, mitad, fondo libre
    Tvals = zeros(numel(Zs), numel(Xs));
    for i = 1:numel(Zs)
       for j = 1:numel(Xs)
           Tvals(i,j) = Tfun(Xs(j), Zs(i));
       end
    end
    disp('Temperaturas [C] (filas Z=1,2,3 m; columnas X=0, e/2, e):');
    disp(Tvals);
    
    % Longitud recomendada por criterio de penetracion Xpen = e
    Zstar = (e/4)^2 * v_max/alpha;
    
    % Flujo de calor promedio y potencia total
    qpp = sqrt(6*alpha*v_med/(pi*L)) * rho*Cp*(Ts-T0);   % [W/m^2]
    QAS = qpp * (b*L);                                   % [W]
    
    fprintf('alpha = %.3e m^2/s\n', alpha);
    fprintf('v_med = %.4f m/s, v_max = %.4f m/s\n', v_med, v_max);
    fprintf('Z* (Xpen=e) = %.4f m\n', Zstar);
    fprintf('q''_AS = %.0f W/m^2,  Q_AS = %.0f W\n', qpp, QAS);

    % grafica
    Zvec  = linspace(0, L, 400);                
    Xpen  = 4*sqrt(alpha*Zvec./v_max);         
    Alt   = L - Zvec;                           
    figure;
    plot(Xpen, Alt, 'r', 'LineWidth', 1.8); hold on;
    xline(e, 'b', 'LineWidth', 1.5);            
    xlabel('Penetración (m)');
    ylabel('Altura (m)');
    title('Penetración térmica a lo largo de la pared');
    xlim([0, 1.05*max(Xpen)]); ylim([0, L]);
    grid off; box on;
%% Ejemplo 5.2.
    function Facf()
    % Datos
    k = 7.0;        % W/mK
    Q = 48.0;       % W
    Tinf = 35.0;    % °C
    L = 0.16;        % m
    D = 0.007;      % m

    % Factores de forma
    Scil = 2*pi*L / log(4*L/D);      % Cilindro semi-infinito
    Sesf = 2*pi*(D/2);               % Esfera semi-infinita
    Spla = 4*(D/2);                  % Placa circular

    % Temperaturas
    Tcil = Tinf + Q/(k*Scil);
    Tesf = Tinf + Q/(k*Sesf);
    Tpla = Tinf + Q/(k*Spla);

    fprintf('Cilindro semi-infinito: T1 = %.2f °C\n', Tcil);
    fprintf('Esfera semi-infinita:  T1 = %.2f °C\n', Tesf);
    fprintf('Placa circular:        T1 = %.2f °C\n', Tpla);

    % Gráfica con variación del diámetro
    Dvals = linspace(0.002, 0.05, 100);
    Tcil_vals = arrayfun(@(d) Tinf + Q/(k*(2*pi*L/log(4*L/d))), Dvals);
    Tesf_vals = arrayfun(@(d) Tinf + Q/(k*(2*pi*(d/2))), Dvals);
    Tpla_vals = arrayfun(@(d) Tinf + Q/(k*(4*(d/2))), Dvals);

    figure;
    plot(Dvals*1000, Tcil_vals, 'B', 'LineWidth', 1.5); hold on;
    plot(Dvals*1000, Tesf_vals, 'k--', 'LineWidth', 1.5);
    plot(Dvals*1000, Tpla_vals, 'c-.', 'LineWidth', 1.5);
    xlabel('Diámetro (mm)');
    ylabel('Temperatura superficial (°C)');
    title('Variación de T con diferentes S');
    legend('Cilindro semi-infinito','Esfera semi-infinita','Placa circular');
    grid on;
end

    %% Ejemplo 5.3.
    clear;
     %Diferencias finitas
    Nx = 4; Ny = 4;   % número de nodos interiores
    Ttol = 1e-4;      % tolerancia
    maxIter = 5000;   % máximo de iteraciones
    % Condiciones de frontera
    Tsup = 100; Tinf = 0; Tizq = 75; Tder = 50;
    T = zeros(Ny+2, Nx+2);
    T(1,:) = Tsup;        
    T(end,:) = Tinf;  
    T(:,1) = Tizq;     
    T(:,end) = Tder; 
    
    for iter = 1:maxIter
       T_viejo = T;
      
       for i = 2:Ny+1
           for j = 2:Nx+1
               T(i,j) = 0.25 * ( T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1) );
           end
       end
       error = max(max(abs(T - T_viejo)));
       if error < Ttol
           fprintf('Convergió en %d iteraciones\n', iter);
           break;
       end
    end
    
    disp('Distribución de temperaturas en nodos interiores:')
    disp(T(2:end-1, 2:end-1))
    figure;
    surf(T);
    title('Distribución de temperaturas');
    xlabel('x'); ylabel('y'); zlabel('T [°C]');
%% Ejemplo 5.4.

    % ----------------------------------------------------------
    % PARTE 1: Parámetros del problema
    % ----------------------------------------------------------
    L = 0.3;         % [m] lado total del bloque
    a = 0.25;        % [m] lado de la parte de acero
    Nx = 101;        % número de nodos en x
    Ny = 101;        % número de nodos en y
    dx = L / (Nx - 1);  % paso en x
    dy = L / (Ny - 1);  % paso en y
    
    % Propiedades térmicas
    k_steel = 15;      % W/m·K  (acero inoxidable)
    k_mat   = 0.25;    % W/m·K  (material aislante)
    
    % Condiciones de contorno y generación
    T_ext   = 40.0;       % °C en el borde exterior
    Q_total = 500.0;      % W generados en el líquido
    V_cav   = a * a * 1.0; % volumen de la cavidad (profundidad unitaria)
    q_gen   = Q_total / V_cav;  % generación volumétrica uniforme [W/m^3]
    
    % ----------------------------------------------------------
    % PARTE 2: Malla y condiciones iniciales
    % ----------------------------------------------------------
    x = linspace(0, L, Nx);
    y = linspace(0, L, Ny);
    [X, Y] = meshgrid(x, y);
    
    % Inicializar temperatura
    T = ones(Ny, Nx) * T_ext;
    
    % Determinar qué nodos pertenecen a la cavidad de acero
    x_cav_min = (L - a) / 2;
    x_cav_max = (L + a) / 2;
    y_cav_min = (L - a) / 2;
    y_cav_max = (L + a) / 2;
    
    % Máscara lógica (1 si el nodo está dentro de la cavidad de acero)
    mask_cav = (X >= x_cav_min) & (X <= x_cav_max) & ...
               (Y >= y_cav_min) & (Y <= y_cav_max);
    
    % ----------------------------------------------------------
    % PARTE 3: Iteracion (Gauss-Seidel)
    % ----------------------------------------------------------
    Ttol = 1e-6;
    max_iter = 10000;
    
    for it = 1:max_iter
        T_old = T;
    
        % Recorremos nodos interiores
        for i = 2:Ny-1
            for j = 2:Nx-1
    
                % Determinar propiedades locales
                if mask_cav(i,j)
                    k = k_steel;
                    q = q_gen;
                else
                    k = k_mat;
                    q = 0.0;
                end
    
                % Ecuación de diferencias finitas con generación
                T(i,j) = ( T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1) ...
                            + q*dx^2/k ) / 4;
            end
        end
    
        % Condiciones de borde (Dirichlet)
        T(1,:)   = T_ext;   % borde superior
        T(end,:) = T_ext;   % borde inferior
        T(:,1)   = T_ext;   % borde izquierdo
        T(:,end) = T_ext;   % borde derecho
    
        % Criterio de convergencia
        err = max(max(abs(T - T_old)));
        if err < Ttol
            fprintf('Convergencia alcanzada en %d iteraciones (error = %.2e)\n', it, err);
            break;
        end
    end
    
    % ----------------------------------------------------------
    % PARTE 4: Resultados numéricos
    % ----------------------------------------------------------
    T_max = max(max(T));
    T_mean_cav = mean(T(mask_cav));
    
    fprintf('q = %.3e W/m^3\n', q_gen);
    fprintf('T_max (dominio) = %.4f °C\n', T_max);
    fprintf('T_media (cavidad) = %.4f °C\n', T_mean_cav);
    
    % ----------------------------------------------------------
    % PARTE 5: Visualización
    % ----------------------------------------------------------
    figure;
    contourf(X, Y, T, 50, 'LineColor', 'none');
    colorbar;
    title('Distribución de temperatura en placa de calentamiento');
    xlabel('x [m]');
    ylabel('y [m]');
    axis equal tight;

    %% Ejemplo 5.5.1
    clc; close all; clear all;
    
    % Conversión inicial en los 5 nodos radiales
    xv = [0 0 0 0 0]';
    
    % Temperatura inicial (K)
    Tv = [873 873 873 873 883]';
    
    % Número de iteraciones
    n = 0;
    
    % Parámetros geométricos
    delz = 0.0945;   % Incremento axial
    M    = 0.179;    
    Mm   = 0.25;     
    
    % Nodos radiales (solo para gráficas)
    j = [1 2 3 4 5];
    
    % Matrices de coeficientes
    
    % Coeficientes para la conversión (ecuación diferencial forma discreta)
    coefix = [ (1-4*Mm)   4*Mm        0           0           0;
               Mm*(1-1/2) (1-2*Mm)   Mm*(1+1/2)   0           0;
               0          Mm*(1-1/4) (1-2*Mm)    Mm*(1+1/4)   0;
               0           0          Mm*(1-1/6) (1-2*Mm)    Mm*(1+1/6);
               0           0          0           Mm*((1-1/8)+(1+1/8)) (1-2*Mm) ];
    
    % Coeficientes para la temperatura
    coefit = [ (1-4*M)    4*M        0           0           0;
               M*(1-1/2)  (1-2*M)   M*(1+1/2)    0           0;
               0          M*(1-1/4) (1-2*M)     M*(1+1/4)    0;
               0           0         M*(1-1/6)  (1-2*M)      M*(1+1/6);
               0           0         0         -0.030        1.030 ];
    
    % Para guardar resultados de cada iteración
    x = zeros(5,11);
    T = zeros(5,11);
    
    while n < 11
    
        % Solución numérica (predicción)
        xnr = coefix * xv;     % conversión predicha
        Tnr = coefit * Tv;     % temperatura predicha
    
        % cálculo de K(T), r_c, r_x, r_T
    
        K = 0.027 .* exp(0.021 .* (Tv - 773));
    
        A = (1.2 .* xv.^2) ./ (K .* (11 + xv).^2);
        B = (1 - xv) ./ (11 + xv);
    
        % Velocidad de reacción (modelo cinético)
        rc = 15100 .* exp(-11000 ./ Tv) .* (B - A);
    
        % Incrementos por reacción
        rx = 169 * delz .* rc;     % cambio en conversión
        rt = -37900 * delz .* rc;  % cambio en temperatura
    
        % Condición especial en nodo 5 (frontera térmica)
        rt(5) = 0;
    
        xn = xnr + rx;
        Tn = Tnr + rt;
    
        xv = xn;
        Tv = Tn;
    
        n = n + 1;
    
        % Guardar resultados
        x(:,n) = xv;
        T(:,n) = Tv;
    
        figure(1)
        hold on
        plot(j, x(:,n), 'o')
        plot(j, x(:,n), 'r')
        grid on
        axis square
        xlabel('Radio (cm)')
        ylabel('Conversión')
    
        figure(2)
        hold on
        plot(j, T(:,n), 'o')
        plot(j, T(:,n), 'r')
        grid on
        axis square
        xlabel('Radio (cm)')
        ylabel('Temperatura (K)')
    end
%% Ejemplo 5.5.2
clc
close all
clear all

% ------------------ Configuración general ------------------

nodosr = [3 5 10 15];        % Número de nodos radiales
col = {'g','b','r','k'};     % Colores para cada malla

% Gráficas
figure(1); hold on; grid on; axis square
xlabel('Radio (cm)')
ylabel('Conversión')
axis([0 5.08 0 0.55])

figure(2); hold on; grid on; axis square
xlabel('Radio (cm)')
ylabel('Temperatura (K)')
axis([0 5.08 800 890])

% ------------------ Bucle sobre mallas radiales ------------------

for t = 1:length(nodosr)

    %-------- Parámetros del reactor --------
    D = 10.16;                  % Diámetro del tubo (cm)
    Nr = nodosr(t);             % Nodos radiales
    Nz = 11;                    % Pasos axiales

    xo  = 0;                    % Conversión inicial
    To  = 873;                  % Temperatura inicial (K)
    To2 = 883;                  % Temperatura en la pared (K)

    delz = 0.0945;              % Paso axial
    delr = D/(2*(Nr-1));        % Paso radial

    M  = 0.179;
    Mm = 0.25;

    % -------- Inicialización --------
    xv = xo*ones(Nr,1);
    Tv = To*ones(Nr,1);
    Tv(Nr) = To2;               % Condición de pared

    coefix = zeros(Nr,Nr);
    coefit = zeros(Nr,Nr);

    % -------- Construcción de matrices de difusión --------
    for i = 1:Nr

        % Nodo central (r = 0)
        if i == 1
            coefix(i,i)   = 1 - 4*Mm;
            coefix(i,i+1) = 4*Mm;

            coefit(i,i)   = 1 - 4*M;
            coefit(i,i+1) = 4*M;
        end

        % Nodos internos
        if i > 1 && i < Nr
            coefix(i,i-1) = Mm*(1 - 1/(2*(i-1)));
            coefix(i,i)   = 1 - 2*Mm;
            coefix(i,i+1) = Mm*(1 + 1/(2*(i-1)));

            coefit(i,i-1) = M*(1 - 1/(2*(i-1)));
            coefit(i,i)   = 1 - 2*M;
            coefit(i,i+1) = M*(1 + 1/(2*(i-1)));
        end

        % Nodo de la pared
        if i == Nr
            coefix(i,i)   = 1 - 2*Mm;
            coefix(i,i-1) = Mm*((1 - 1/(2*(i-1))) + (1 + 1/(2*(i-1))));

            coefit(i,i)   = 1.030;
            coefit(i,i-1) = -0.030;
        end
    end

    % Coordenada radial
    r = (0:Nr-1)*delr;

    % -------- Bucle axial --------
    for n = 1:Nz

        % Difusión radial
        xnr = coefix * xv;
        Tnr = coefit * Tv;

        % Cinética
        K  = 0.027 * exp(0.021*(Tv - 773));
        A  = (1.2*xv.^2) ./ (K .* (11 + xv).^2);
        B  = (1 - xv) ./ (11 + xv);
        rc = 15100 * exp(-11000 ./ Tv) .* (B - A);

        % Términos de reacción
        rx =  169 * delz * rc;
        rt = -37900 * delz * rc;
        rt(Nr) = 0;              % Condición térmica en la pared

        % Actualización
        xv = xnr + rx;
        Tv = Tnr + rt;

        % Gráficas
        figure(1)
        plot(r, xv, col{t})

        figure(2)
        plot(r, Tv, col{t})

    end
end
%% Ejemplo 5.5.3

clc
close all
clear all

% PARÁMETROS GENERALES

nodosr = [5];                          % Nodos en el radio
col = ['y','g','b','r','k'];           % Colores para grafica

% MODELO POR DIFERENCIAS FINITAS
for t = 1:length(nodosr)
  
    % GEOMETRÍA Y MALLA
    D  = 10.16;                        % Diámetro del tubo [cm]
    r  = nodosr(t);                    % Nodos radiales
    s  = 11;                           % Nodos axiales
    xo = 0;                            % Conversión inicial
    To = 873;                          % Temperatura inicial [K]
    To2 = 893;                         % Temperatura en la pared [K]
    
    delz = 0.0945;                     % Paso axial
    delr = D/(2*(r-1));                % Paso radial
    
    % PARÁMETROS ADIMENSIONALES
    M  = 0.179;                        % Energía
    Mm = 0.25;                         % Masa
    
    % CONDICIONES INICIALES
    x(:,1) = [xo xo xo xo xo]';         % Conversión
    T(:,1) = [To To To To To2]';        % Temperatura
    
    % INICIALIZACIÓN
    coefix = zeros(r,r);               % Matriz difusión masa
    coefit = zeros(r,r);               % Matriz difusión energía
    
    for i = 1:r
        xv(i,1) = xo;
        Tv(i,1) = To;
        
        for j = 1:r
            k = i;
            p = j;
            
            % NODO CENTRAL (r = 0)
            if k == 1
                if p == 1
                    coefix(i,j)   = (1 - 4*Mm);
                    coefit(i,j)   = (1 - 4*M);
                    coefix(i,j+1) = 4*Mm;
                    coefit(i,j+1) = 4*M;
                end
            end
            
            % NODOS INTERNOS
            if k > 1
                if k < r
                    if k == p
                        coefix(i,j)   = (1 - 2*Mm);
                        coefix(i,j-1) = Mm*(1 - 1/(2*(i-1)));
                        coefix(i,j+1) = Mm*(1 + 1/(2*(i-1)));
                        
                        coefit(i,j)   = (1 - 2*M);
                        coefit(i,j-1) = M*(1 - 1/(2*(i-1)));
                        coefit(i,j+1) = M*(1 + 1/(2*(i-1)));
                    end
                end
            end
            
            % NODO DE PARED
            if k == r
                Tv(i,1) = To2;
                if p == r
                    coefix(i,j)   = (1 - 2*Mm);
                    coefix(i,j-1) = Mm*((1 - 1/(2*(i-1))) + (1 + 1/(2*(i-1))));
                    
                    coefit(i,j)   = 1.030;
                    coefit(i,j-1) = -0.030;
                end
            end
        end
    end
    
    % INTEGRACIÓN AXIAL
    n = 1;
    while n < s
        
        % Difusión radial
        xnr = coefix * xv;
        Tnr = coefit * Tv;
        
        % Cinética
        K = 0.027 * exp(0.021*(Tv - 773));
        A = (1.2*xv.^2) ./ (K.*(11 + xv).^2);
        B = (1 - xv) ./ (11 + xv);
        rc = 15100 * exp(-11000./Tv) .* (B - A);
        
        % Términos fuente
        rx = 169 * delz * rc;
        rt = -37900 * delz * rc;
        rt(r) = 0;                      % Condición de pared
        
        % Actualización
        xn = xnr + rx;
        Tn = Tnr + rt;
        
        xv = xn;
        Tv = Tn;
        
        n = n + 1;
        x(:,n) = xn;
        T(:,n) = Tn;
    end
    
    % EJE AXIAL
    j = [0:(s-1)] * delz;
    
    % GRÁFICAS
    for h = 1:r
        figure(1)
        hold on
        plot(j, x(h,:), 'o')
        plot(j, x(h,:), col(h))
        grid on
        
        figure(2)
        hold on
        plot(j, T(h,:), 'o')
        plot(j, T(h,:), col(h))
        grid on
    end
    
    % RESULTADOS FINALES
    xfin = x
    Tfin = T
    nodosr = [5];
end

% MODELO CON ODE45
nodosr = [5];
col = ['y','g','b','r','k'];

for t = 1:length(nodosr)
    
    D  = 0.1016;                       % Diámetro [m]
    r  = nodosr(t);
    s  = 11;
    xo = 0;
    To = 873;
    
    delz = 0.0945;
    L = s * delz;
    delr = D/(2*(r-1));
    
    M  = 0.179;
    Mm = 0.25;
    
    [z,x] = ode45(@bidimensional,[0 L], ...
        [xo xo xo xo xo To To To To To+20],[],Mm,M,delr,delz);
    
    figure(1)
    hold on
    plot(z,x(:,1),col(1))
    plot(z,x(:,2),col(2))
    plot(z,x(:,3),col(3))
    plot(z,x(:,4),col(4))
    plot(z,x(:,5),col(5))
    xlabel('Longitud')
    ylabel('Conversion cada nodo')
    grid on
    axis square
    
    figure(2)
    hold on
    plot(z,x(:,6),col(1))
    plot(z,x(:,7),col(2))
    plot(z,x(:,8),col(3))
    plot(z,x(:,9),col(4))
    plot(z,x(:,10),col(5))
    xlabel('Longitud')
    ylabel('Temperatura cada nodo (K)')
    grid on
    axis square
end

%% Funciones 

    function dYdz = bidimensional(z, Y, Mm, M, delr, delz)
    % Y es un vector que contiene [x1, x2, x3, x4, x5, T1, T2, T3, T4, T5]
    % donde 1-5 son nodos radiales de conversion y 6-10 son de temperatura
    
    dYdz = zeros(10,1);
    r = 5; % Número de nodos
    
    % Extraer vectores para facilitar lectura
    xv = Y(1:r);
    Tv = Y(r+1:end);
    
    % --- Cálculo de la Cinética ---
    K = 0.027 * exp(0.021 * (Tv - 773));
    A = (1.2 * xv.^2) ./ (K .* (11 + xv).^2);
    B = (1 - xv) ./ (11 + xv);
    rc = 15100 * exp(-11000 ./ Tv) .* (B - A);
    
    % Tasas de cambio (reacción)
    rx = 169 * rc;
    rt = -37900 * rc;
    rt(r) = 0; % Condición en la pared para la temperatura
    
    % --- Diferenciación Espacial (Balance de Masa y Energía) ---
    for i = 1:r
        % Coeficientes de difusión/conducción radial simplificados para ODE
        if i == 1 % Centro del tubo (Simetría)
            dif_x = 4 * Mm * (xv(2) - xv(1));
            dif_t = 4 * M * (Tv(2) - Tv(1));
        elseif i < r % Nodos intermedios
            dif_x = Mm * ((1 - 1/(2*(i-1)))*xv(i-1) - 2*xv(i) + (1 + 1/(2*(i-1)))*xv(i+1));
            dif_t = M * ((1 - 1/(2*(i-1)))*Tv(i-1) - 2*Tv(i) + (1 + 1/(2*(i-1)))*Tv(i+1));
        else % Pared del tubo
            dif_x = Mm * ((2)*xv(i-1) - 2*xv(i));
            dif_t = 0.030 * (Tv(i-1) - Tv(i)); % Basado en tus coeficientes del clear anterior
        end
        
        % dX/dz y dT/dz
        dYdz(i) = (dif_x / delz) + rx(i);
        dYdz(i+r) = (dif_t / delz) + rt(i);
    end
end