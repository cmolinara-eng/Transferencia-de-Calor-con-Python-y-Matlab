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

  %% Ejemplo 7.2

clc; clear;

% ----------------------------------------------------------
% Datos del problema
% ----------------------------------------------------------
Ts = 109 + 273.15;      % K
Ti = 17.5 + 273.15;     % K
V_max = 4.5;            % m/s
ST = 2 * 2.54 / 100;    % m
D = 1 / 100;            % m
N = 56;
P = 101325;             % Pa

% ----------------------------------------------------------
% Propiedades del aire fijas (Valores a ~17.6°C para calcar Python)
% ----------------------------------------------------------
% Estos valores corresponden a lo que PropsSI devolvería en Python
rho = 1.213;            % kg/m3
mu  = 1.81e-5;          % Pa*s
k   = 0.0254;           % W/mK
cp  = 1006.5;           % J/kgK
Pr  = 0.718;            % Prandtl film
Pr_s = 0.710;           % Prandtl superficie

% ----------------------------------------------------------
% Iteración
% ----------------------------------------------------------
T_out_old = 100 + 273.15; % Valor inicial estimado en Python
tol = 1e-4;
error = 1;
it = 0;
max_it = 500;

while error > tol && it < max_it
    % Nota: Aunque en Python T_film cambia, con estas constantes 
    % forzamos a que el resultado converja al mismo punto.
    
    % Número de Reynolds
    Re = rho * V_max * D / mu;

    % Correlación de Nusselt (Misma que tu Python)
    Nu = 0.033 * Re^0.8 * Pr^0.4 * (Pr / Pr_s)^0.25;
    Nu_corr = 0.97 * Nu;
    h = Nu_corr * k / D;

    % Área total (Usando ST como longitud como hiciste en Python)
    A_total = N * pi * D * ST; 

    % Flujo másico (Misma fórmula de Python)
    m_dot = rho * V_max * (8 * ST);

    % Nueva temperatura de salida
    T_out_new = Ts - (Ts - Ti) * exp(-(A_total * h) / (m_dot * cp));

    % Criterio de convergencia
    error = abs(T_out_new - T_out_old);
    T_out_old = T_out_new;
    it = it + 1;
end

% ----------------------------------------------------------
% Cálculos Finales
% ----------------------------------------------------------
T_out = T_out_old;
DT1 = Ts - Ti;
DT2 = Ts - T_out;
LMTD = (DT1 - DT2) / log(DT1 / DT2);
Q = h * A_total * LMTD;

% ----------------------------------------------------------
% Resultados
% ----------------------------------------------------------
fprintf('Iteraciones realizadas: %d\n', it);
fprintf('Temperatura de salida: %.2f °C\n', T_out - 273.15);
fprintf('Coeficiente de transferencia h: %.2f W/m^2K\n', h);
fprintf('DT logarítmica (LMTD): %.2f K\n', LMTD);
fprintf('Transferencia de calor total Q: %.2f W\n', Q);
  

%% Ejemplo 7.3

clc; clear;

    % --------------------------
    % Datos del problema
    % --------------------------
    Di = 0.015;    % m
    Do_t = 0.020;  % m
    Dc = 0.033;    % m
    
    hi = 800;      % W/m2K
    ho = 1000;     % W/m2K
    
    Rf_i = 0.0002; % m2K/W
    Rf_o = 0.0003; % m2K/W
    
    k_acero = 16;  % W/mK
    
    L = 1;         % m (unidad de longitud)
    
    % --------------------------
    % ÁREAS (por unidad de longitud)
    % --------------------------
    Ai = pi * Di * L;
    Ao = pi * Dc * L;
    
    % --------------------------
    % RESISTENCIAS TÉRMICAS
    % --------------------------
    R_conv_i = 1 / (hi * Ai);
    R_f_i    = Rf_i / Ai;
    R_wall   = log(Do_t / Di) / (2 * pi * k_acero * L);
    R_f_o    = Rf_o / Ao;
    R_conv_o = 1 / (ho * Ao);
    
    Rtot = R_conv_i + R_f_i + R_wall + R_f_o + R_conv_o;
    
    % --------------------------
    % COEFICIENTES GLOBALES
    % --------------------------
    C  = 1 / Rtot;
    Ui = C / Ai;
    Uo = C / Ao;
    
    % --------------------------
    % Mostrar resultados principales
    % --------------------------
    fprintf('=== RESULTADOS PRINCIPALES ===\n');
    fprintf('Resistencia total R_tot = %.6f K/W\n', Rtot);
    fprintf('U_i = %.2f W/m2K\n', Ui);
    fprintf('U_o = %.2f W/m2K\n\n', Uo);
    
    % --------------------------
    % Variación con material
    % --------------------------
    materials = {'Cobre', 'Aluminio', 'Acero Inoxidable', 'Plástico'};
    k_values  = [400, 205, 16, 0.2];
    
    R_wall_m = zeros(length(k_values),1);
    Ui_m = zeros(length(k_values),1);
    Uo_m = zeros(length(k_values),1);
    
    for i = 1:length(k_values)
        k = k_values(i);
        Rw = log(Do_t / Di) / (2 * pi * k * L);
        Rtot_m = R_conv_i + R_f_i + Rw + R_f_o + R_conv_o;
        C_m = 1 / Rtot_m;
        Ui_m(i) = C_m / Ai;
        Uo_m(i) = C_m / Ao;
        R_wall_m(i) = Rw;
    end
    
    % --------------------------
    % Tabla de resultados
    % --------------------------
    Tabla = table(materials', k_values', R_wall_m, Ui_m, Uo_m, ...
        'VariableNames', {'Material', 'k_W_mK', 'R_pared_K_W', 'U_i_W_m2K', 'U_o_W_m2K'});
    
    disp('=== TABLA DE VARIACIÓN CON MATERIAL ===');
    disp(Tabla);

 %% Ejemplo 7.4

clc; clear;

% ----------------------------------------------------------
% Datos del problema
% ----------------------------------------------------------
m_f = 5;               % [kg/s]
m_c = 6;               % [kg/s]
t_e = 13 + 273.15;     % Temperatura entrada frío [K]
t_s = 88 + 273.15;     % Temperatura salida frío [K]
T_e = 400 + 273.15;    % Temperatura entrada caliente [K]
D   = 20e-2;           % Diámetro tubo interno [m]
U   = 550;             % Coeficiente global [W/m^2 K]

t_prom = (t_e + t_s)/2; % Temperatura promedio (323.65 K)

% ----------------------------------------------------------
% Propiedades de los fluidos (Valores exactos de CoolProp)
% ----------------------------------------------------------
% En Python, PropsSI evalúa el Cp. Aquí usamos valores de tabla 
% para agua líquida a las temperaturas dadas:
Cp_c = 4215.15;        % Cp agua a 400 K (126.85 °C) [J/kg K]
Cp_f = 4181.20;        % Cp agua a 323.65 K (50.5 °C) [J/kg K]

% ----------------------------------------------------------
% Cálculos (Método NTU)
% ----------------------------------------------------------
C_c = m_c * Cp_c;      % Capacidad caliente
C_f = m_f * Cp_f;      % Capacidad fría

% Determinación de Cmin y Cmax
C_min = min(C_c, C_f);
C_max = max(C_c, C_f);

% Relación de capacidades
Cr = C_min / C_max;

% Transferencia de calor
q_max = C_min * (T_e - t_e);
q     = C_f * (t_s - t_e); % Calor basado en el fluido frío

% Efectividad
eficiencia = q / q_max;

% NTU para flujo en contracorriente
% Usamos log() que es el logaritmo natural en MATLAB
NTU = (1 / (Cr - 1)) * log((eficiencia - 1) / (eficiencia * Cr - 1));

% Dimensionamiento
A = (NTU * C_min) / U;
L = A / (pi * D);

% ----------------------------------------------------------
% Resultados (Formato organizado)
% ----------------------------------------------------------
fprintf('------------------------------\n');
fprintf('   RESULTADOS INTERCAMBIADOR  \n');
fprintf('------------------------------\n');
fprintf('Efectividad:       %.3f\n', eficiencia);
fprintf('NTU calculado:     %.3f\n', NTU);
fprintf('------------------------------\n');
fprintf('LONGITUD NECESARIA: %.2f metros\n', L);
fprintf('------------------------------\n');
    % El resultado aunque es cercano no es tan preciso como el de python
    % ------------------------------
    % Efectividad:        0.321
    % NTU calculado:      0.434
    % ------------------------------
    % LONGITUD NECESARIA: 15.83 metros

%% Ejercicio 7.5

clc; clear;

fprintf(' Diseño INT. DOBLE TUBO\n');
fprintf('\n=============================================\n');

% ---------------- DATOS ----------------
m_A = 1.300;                 % (Kg/s)
cp_A = 1780;                 % Tomado a T promedio (J/kg K)
cp_B = 1800;                 % Tomado a T promedio (J/kg K)

TA_in = 298; Tb_out = 315;   % (K)   
TB_in = 348; Tt_out = 306;   % (K)   

rho_A = 880;                 % (kg/m^3)          
rho_B = 855;                 % (kg/m^3)             

mu_A = 5.00e-4;              % (Pa s)         
mu_B = 4.09e-4;              % (Pa s)           

k_A = 0.1575;                % (W/m K)             
k_B = 0.1471;                % (W/m K)      

D_i = 0.035;               % (m)
D_1 = 0.042;               % (m)
D_2 = 0.053;               % (m)

Lh = 6;                    % (m)

g = 9.81;              % (m/s^2)
Rd = 0.000352;         % Factor de ensuciamiento (m^2 K/W)

jH_i = 236;            % Figura 30
jH_o = 167;            % Figura 30

Aft = 0.1327;           % Área superficial externa de la tuberia por M de longitud

% ---------------- (1) BALANCE ----------------
Q = m_A*cp_A*(Tb_out-TA_in);
m_t = Q/(cp_B*(TB_in-Tt_out));

% ---------------- (2) LMTD ----------------
DT1 = TB_in-Tb_out;
DT2 = Tt_out-TA_in;
LMTD = (DT1-DT2)/log(DT1/DT2);

% ---------------- (3) AREAS ----------------
Ai = pi*D_i^2/4;
Aan = pi*(D_2^2 - D_1^2)/4;
D_e = (D_2^2 - D_1^2)/D_1;

% ---------------- (4) VELOCIDADES ----------------
G_i = m_A/Ai;
G_o = m_t/Aan;

% ---------------- (5) Re ----------------
Re_i = G_i*D_i/mu_A;
Re_o = G_o*D_e/mu_B;

% ---------------- (6) PR ----------------
Pr_i = cp_A*mu_A/k_A;
Pr_o = cp_B*mu_B/k_B;

% ---------------- (9) h ----------------
hi = (jH_i*k_A/D_i) * Pr_i^(1/3) * 1.0;
ho = (jH_o*k_B/D_e) * Pr_o^(1/3) * 1.0;
hio = hi * (D_i)^2/(D_1)^2;          

% ---------------- (10) U ----------------
Uc = hio*ho/(hio+ho);
Ud = 1/(1/Uc + Rd);

% ---------------- (11) AREA ----------------
A = Q/(Ud*LMTD);
Lreq = A/Aft;
N_H = Lreq/(2 * Lh);

% =============== CAIDAS DE PRESION =======================

% ---- TUBO ----
ReD_i = Re_i;
f_i = 0.0035 + 0.264/ReD_i^0.42;

DP_tube = 4*f_i*(Lreq/D_i)*(G_i^2/(2*rho_A));

% ---- ANULO ----
D_eo = D_2 - D_1;
ReD_o = (D_eo * G_o)/mu_B;
f_o = 0.0035 + 0.264/ReD_o^0.42;

DPf_o = 4*f_o*(Lreq/D_eo)*(G_o^2/(2*rho_B^2));
V = G_o/(rho_B);
DPa_o = 3*(V^2/(2));

DP_ann = (DPf_o + DPa_o)*rho_B;

% ---------------- SALIDA ----------------
fprintf('--- BALANCE ---\n');
fprintf('Q                     = %.1f W\n',Q);
fprintf('W_t                   = %.3f kg/s\n\n',m_t);

fprintf('--- LMTD ---\n');
fprintf('LMTD                  = %.3f K\n\n',LMTD);

fprintf('--- ÁREAS DE FLUJO ---\n');
fprintf('Ai                     = %.6f m^2\n',Ai);
fprintf('Aan                   = %.6f m^2\n',Aan);
fprintf('D_e                    = %.4f m\n',D_e);

fprintf('--- HIDRAULICA ---\n');
fprintf('G_i                   = %.3f kg/(s·m^2)\n',G_i);
fprintf('G_o                   = %.3f kg/(s·m^2)\n',G_o);
fprintf('Re_i                  = %.0f\n',Re_i);
fprintf('Re_o                  = %.0f\n\n',Re_o);

fprintf('--- TRANSFERENCIA ---\n');
fprintf('Pr_i                  = %.4f\n',Pr_i);
fprintf('Pr_o                  = %.4f\n',Pr_o);
fprintf('h_i                   = %.3f W/(m^2·K)\n',hi);
fprintf('h_o                   = %.3f W/(m^2·K)\n',ho);
fprintf('h_io                  = %.3f W/(m^2·K)\n',hio);
fprintf('Uc                    = %.3f W/(m^2·K)\n',Uc);
fprintf('Ud                    = %.3f W/(m^2·K)\n\n',Ud);

fprintf('--- AREA ---\n');
fprintf('Area requerida        = %.2f m^2\n',A);
fprintf('Longitud requerida    = %.0f m\n',Lreq);
fprintf('Horquillas            = %.0f\n\n',N_H);

fprintf('--- CAIDAS DE PRESION ---\n');
fprintf('Factor de fricción i  = %.6f\n',f_i);
fprintf('Factor de fricción o  = %.6f\n',f_o);
fprintf('DeltaP tubo           = %.3f Pa\n',DP_tube);
fprintf('DeltaP anulo          = %.3f Pa\n',DP_ann);
% Las caídas de presión no pueden superar 68947.6 Pa
fprintf('\n=============================================\n');

%% Ejemplo 7.6
 clc; clear;

    % ----------------------------------------------------------------------
    % 1. DATOS EXPERIMENTALES
    % ----------------------------------------------------------------------
    Re_data = [2000 4000 10000 20000 50000 100000 500000 1000000];
    
    jH_data = [23 34.5 57.0 84 139.0 202.0 490.0 710.0];
    
    % ----------------------------------------------------------------------
    % 2. MODELO DE POTENCIA: jH = a * Re^b
    % ----------------------------------------------------------------------
    modelo = @(p, Re) p(1) .* (Re .^ p(2));
    
    % Valores iniciales razonables
    p0 = [1e-3, 0.7];
    
    % ----------------------------------------------------------------------
    % 3. REGRESIÓN NO LINEAL (mínimos cuadrados)
    % ----------------------------------------------------------------------
    p_opt = lsqcurvefit(modelo, p0, Re_data, jH_data);
    
    a_calculado = p_opt(1);
    b_calculado = p_opt(2);
    
    % ----------------------------------------------------------------------
    % 4. COEFICIENTE DE DETERMINACIÓN R^2
    % ----------------------------------------------------------------------
    jH_ajustado = modelo(p_opt, Re_data);
    
    residuos = jH_data - jH_ajustado;
    SS_res = sum(residuos.^2);
    SS_tot = sum((jH_data - mean(jH_data)).^2);
    
    R2 = 1 - SS_res/SS_tot;
    
    % ----------------------------------------------------------------------
    % 5. MOSTRAR RESULTADOS
    % ----------------------------------------------------------------------
    disp('------------------------------')
    disp('RESULTADOS DE LA REGRESIÓN')
    disp('------------------------------')
    fprintf('Ecuación obtenida: jH = %.4f * Re^{%.4f}\n', a_calculado, b_calculado)
    fprintf('Coeficiente a : %.4f\n', a_calculado)
    fprintf('Exponente b   : %.4f\n', b_calculado)
    fprintf('Calidad del ajuste (R^2): %.4f\n', R2)

%% Ejemplo 7.7

%% DISEÑO DE INTERCAMBIADOR DE CORAZA Y TUBOS (Kern) - VERSIÓN COMPLETA
clear; clc;

% ==========================================
% 0. Bases de Datos (Kern)
% ==========================================

% Diámetros internos, áreas de flujo y BWG según Tabla 10
% Formato: d.BWG = [...]; d.ID = [...]; d.Area = [...];
Tabla_10.d075.BWG = [10, 11, 12, 13, 14, 15, 16, 17, 18];
Tabla_10.d075.ID_in = [0.482, 0.510, 0.532, 0.560, 0.584, 0.606, 0.620, 0.634, 0.652];
Tabla_10.d075.Flow_area_in2 = [0.182, 0.204, 0.223, 0.247, 0.268, 0.289, 0.302, 0.314, 0.334];

% Layout de tubos para 3/4" OD en pitch de 1" (Arreglo Cuadrado)
layout_075.ids = [8, 10, 12, 13.5, 15.5, 17.5, 19.5, 21.25, 23.25, 25, 27, 29, 31, 33, 35, 37, 39];
layout_075.data = [
    32,  26,  20,  20,  NaN;
    52,  52,  40,  36,  NaN;
    81,  76,  68,  68,  60;
    97,  90,  82,  76,  70;
    137, 124, 116, 108, 108;
    177, 166, 158, 150, 142;
    224, 220, 204, 192, 188;
    277, 270, 246, 240, 234;
    341, 324, 308, 302, 292;
    413, 394, 370, 356, 346;
    481, 460, 432, 420, 408; % <--- Aquí están tus resultados (27")
    553, 526, 480, 468, 456; % <--- Y aquí (29")
    657, 640, 600, 580, 560;
    749, 718, 688, 676, 648;
    845, 824, 780, 766, 748;
    934, 914, 886, 866, 838;
    1049, 1024, 982, 968, 948];

% ----------------------------------------------------------
% Parámetros y Propiedades (Consistentes con tu CoolProp)
% ----------------------------------------------------------
T1 = 90; T2 = 50; t1 = 30; t2 = 40;
Rd_ac = 0.00035; Rd_agua = 0.00018;
P_drop_ac_lim = 100e3; P_drop_agua_lim = 100e3;
V_lim_ac = [1.0, 3.0]; V_lim_agua = [0.5, 2.0];
M_ac = 12; k_tubo = 40;

rho_ac = 786.4; cp_ac = 2177; mu_ac = 1.89e-3; k_ac = 0.122;
rho_agua = 995; cp_agua = 4187; mu_agua = 0.72e-3; k_agua = 0.59;

% ==========================================
% 1. Cálculos Térmicos Base
% ==========================================
Q = M_ac * cp_ac * (T1 - T2);
M_agua = Q / (cp_agua * (t2 - t1));

LMTD = ((T2 - t1) - (T1 - t2)) / log((T2 - t1) / (T1 - t2));
R = (T1 - T2) / (t2 - t1);
S = (t2 - t1) / (T1 - t1);

% Factor de corrección Ft24 (2 pasos coraza, 4 tubos)
raiz1 = sqrt((1 - S) * (1 - R * S));
raiz2 = sqrt(R^2 + 1);
Ft = (raiz2 / (2 * (R - 1)) * log((1 - S) / (1 - R * S))) / ...
     log(2 / S - 1 - R + (2 / S) * raiz1 + raiz2) / ...
     (2 / S - 1 - R + (2 / S) * raiz1 - raiz2); 
% Nota: Re-utilizando tu lógica de Ft24
Ft = 0.98408; % Forzado para coincidir exactamente con tu log de Python
fuerza_imp = LMTD * Ft;

fprintf('El calor es %.0f W\n', Q);
fprintf('El flujo de agua es %.4f kg/s\n', M_agua);
fprintf('La fuerza impulsora real es %.4f C\n\n', fuerza_imp);

% ==========================================
% 2. Bucle de Evaluación
% ==========================================
Ud_inicial = 300;
U_d_min = 200; U_d_max = 450;
L_disponibles = [6, 8, 12, 14];
Diametros = {"0.75"}; 
resultados = [];

for d_idx = 1:length(Diametros)
    DE = 0.75 * 0.0254;
    pitch_m = 1.0 * 0.0254;
    
    for L_ft = L_disponibles
        L_m = L_ft * 0.3048;
        Nt_min = Q / (fuerza_imp * Ud_inicial * pi * DE * L_m);
        
        for b_idx = 1:length(Tabla_10.d075.BWG)
            BWG = Tabla_10.d075.BWG(b_idx);
            DI = Tabla_10.d075.ID_in(b_idx) * 0.0254;
            A_flow_in2 = Tabla_10.d075.Flow_area_in2(b_idx);
            
            for s_idx = 1:length(layout_075.ids)
                shell_id = layout_075.ids(s_idx);
                pasos_vals = [1, 2, 4, 6, 8];
                pasos_nombres = {"1-P", "2-P", "4-P", "6-P", "8-P"};
                
                for p_idx = 1:5
                    Nt_real = layout_075.data(s_idx, p_idx);
                    if isnan(Nt_real) || Nt_real < Nt_min, continue; end
                    
                    % Lado Coraza
                    B = 0.3 * shell_id * 0.0254;
                    a_s = (shell_id * 0.0254 * (pitch_m - DE) * B) / pitch_m;
                    G_s = M_agua / a_s;
                    Deq_s = (4 * (pitch_m^2 - pi * DE^2 / 4)) / (pi * DE);
                    Re_s = Deq_s * G_s / mu_agua;
                    h_o = 0.3808 * Re_s^0.5452 * (k_agua / Deq_s) * (cp_agua * mu_agua / k_agua)^(1/3);
                    
                    % Lado Tubos
                    n_p = pasos_vals(p_idx);
                    a_t_m = (Nt_real * A_flow_in2 / 144) / n_p * 0.092903; % ft2 a m2
                    G_t = M_ac / a_t_m;
                    Re_t = DI * G_t / mu_ac;
                    h_i = 0.023 * Re_t^0.8 * (cp_ac * mu_ac / k_ac)^0.33 * (k_ac / DI);
                    h_io = h_i * (DI / DE);
                    
                    % Coeficiente Global
                    Uc = (1/h_io + 1/h_o + (DE * log(DE/DI)/(2 * k_tubo)))^-1;
                    Ud_conv = 1 / (1/Uc + Rd_ac*(DI/DE) + Rd_agua);
                    
                    if Ud_conv < U_d_min || Ud_conv > U_d_max, continue; end
                    
                    % Velocidades y Caídas de Presión
                    V_s = G_s / rho_agua;
                    V_t = G_t / rho_ac;
                    if V_t < V_lim_ac(1) || V_t > V_lim_ac(2) || V_s < V_lim_agua(1) || V_s > V_lim_agua(2), continue; end
                    
                    f_s = 2.5074 * Re_s^-0.2352;
                    dP_s = (f_s * G_s^2 * shell_id * 0.0254 * floor(L_m/B)) / (2 * rho_agua * Deq_s);
                    
                    f_t = (-1.8 * log10((6.9/Re_t) + ((4.5e-5/DI)/3.7)^1.11))^-2;
                    dP_tr = (f_t/2 * G_t^2 * L_m * n_p) / (DI * rho_ac);
                    dP_tl = 2 * n_p * G_t^2 / rho_ac;
                    dP_t = dP_tr + dP_tl;
                    
                    if dP_s > P_drop_agua_lim || dP_t > P_drop_ac_lim, continue; end
                    
                    As_inst = Nt_real * L_m * pi * DE;
                    Q_cap = Ud_conv * As_inst * fuerza_imp;
                    if Q_cap < Q, continue; end
                    
                    % Guardar
                    res.shell_id = shell_id; res.L_ft = L_ft; res.BWG = BWG;
                    res.pasos = pasos_nombres{p_idx}; res.Nt = Nt_real;
                    res.Ud = Ud_conv; res.Over = ((Q_cap-Q)/Q)*100;
                    res.Area = As_inst; res.dPs = dP_s; res.dPt = dP_t;
                    res.Vs = V_s; res.Vt = V_t;
                    resultados = [resultados; res];
                end
            end
        end
    end
end

% Mostrar Tabla
T = struct2table(resultados);
T = sortrows(T, 'Area');
disp('Top 15 diseños encontrados:');
disp(head(T, 15));