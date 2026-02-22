%%Ejemplo 6.1

function evolucion_temperatura_varilla()
clc; clear; close all;
    
k = 15;
rho = 7000;
cp = 500;
alpha = k / (rho * cp);
L = 0.4;
T1 = 200;
Ts = 250;
T2 = 123;
    
N_terms = 250;
x = linspace(0, L, 600);
t = linspace(0.1, 5000, 700);
[X, T] = meshgrid(x, t);
Z = zeros(size(X));
    
for i = 1:length(t)
    t_val = t(i);
    Z(i,:) = T_xt_matlab(x, t_val, N_terms, L, T1, Ts, T2, alpha);
end
    
figure('Name', 'Evolución de T(x,t) en la varilla', 'NumberTitle', 'off', 'Position', [100 100 800 600]);
surf(X, T, Z, 'EdgeColor', 'none', 'FaceColor', 'interp'); 
    
xlabel('Posición x [m]');
ylabel('Tiempo t [s]');
zlabel('Temperatura [K]');
title('Evolución de T(x,t) en la varilla');
view(45, 30);
colormap('parula'); 
grid on;
axis tight;
colorbar;
end
    
function T_x = T_xt_matlab(x, t, N_terms, L, T1, Ts, T2, alpha)
    
    theta = zeros(size(x));
        
    for n = 1:N_terms
        term1 = (4/pi) * (Ts - T2) * ((-1)^(n-1)) / (2*n - 1);
        term2 = (8 * (Ts - T1)) / ((2*n - 1)^2 * pi^2);
        An = term1 - term2;
            
        lam = (2*n - 1) * pi / (2*L);
            
        theta = theta + An * cos(lam * x) * exp(-(lam^2 * alpha * t));
    end
        
    T_x = theta + T2;
end

%% Ejempo 6.2

clc; clear; close all
    
Ts=400; T0=30; k=50; Den=7980; Cp=510; R=0.02;
dif=k/(Den*Cp);
    
% Estimación de convergencia en el centro
r=0; tiempo=0; n=5000; f=1:n; Tet=0;
for j=1:n
    Lam1=(4*j-1)*(pi/4)/R;
    dTet=(2/R)*(besselj(0,Lam1*r)*exp(-dif*Lam1^2*tiempo))...
                /(Lam1*besselj(1,Lam1*R));
    Tet=Tet+dTet;
    Val(j)=Tet;
end
Temp=Ts+Tet*(T0-Ts);
    
figure(1)
plot(f,Val)
xlabel('Términos de sumatoria'); ylabel('Teta')
title('Convergencia de la serie en el centro')
% Superficie de temperaturas
n=200; nt=60; nr=50;
t=linspace(0,nt,80); r=linspace(0,R,nr);
[Rm,T]=meshgrid(r,t); Temp=zeros(size(Rm));
    
for it=1:length(t)
    for ir=1:length(r)
        Tet=0;
        for j=1:n
            Lam1=(4*j-1)*(pi/4)/R;
            dTet=(2/R)*(besselj(0,Lam1*r(ir))*exp(-dif*Lam1^2*t(it)))...
                          /(Lam1*besselj(1,Lam1*R));
            Tet=Tet+dTet;
        end
        Temp(it,ir)=Ts+Tet*(T0-Ts);
    end
end
    
figure(2)
surf(Rm,T,Temp,'FaceAlpha',0.6,'EdgeColor','interp') 
xlabel('Radio (m)')
ylabel('Tiempo (s)')
zlabel('Temperatura (°C)')
title('Perfiles de temperatura radiales con respecto al tiempo (t=60 s)')
grid on
colormap(jet)

%% Ejemplo 6.3
% Enfriamiento de esfera
clear; close all; clc;
% Datos del problema
d = 0.02;            % m
k = 195;             % W/(m K)
h = 30;              % W/(m^2 K)
rho = 2000;          % kg/m^3
cp = 0.9e3;        % J/(kg K)  (0.9 kJ/kgK)
T0 = 300;            % °C
Tinf = 90;          % °C
Tt = 130;       % °C objetivo
% Geometría y propiedades derivadas
Lc = d/6;                    % m (V/A para esfera)
V = pi * d^3 / 6;
A = pi * d^2;
m = rho * V;
% Cálculos
Bi = h * Lc / k;
tau = m * cp / (h * A);      % s
t_req = - tau * log( (Tt - Tinf) / (T0 - Tinf) );  % s
% Mostrar resultados
fprintf(' Bi = %.4e\\n', Bi);
fprintf(' tau = %.2f s (%.2f min)\\n', tau, tau/60);
fprintf(' Tiempo para alcanzar %.1f °C: t = %.2f s (%.2f min)\\n', Tt, t_req, t_req/60);
% Graficar la evolución de la temperatura y marcar el tiempo objetivo
tmax = max(1.5*t_req, 1200); % s, ventana suficiente
tvec = linspace(0, tmax, 500);
Tvec = Tinf + (T0 - Tinf) * exp(-tvec / tau);
figure('Units','normalized','Position',[0.1 0.15 0.6 0.5])
plot(tvec, Tvec, 'LineWidth', 1.6)
hold on
yline(Tt, '--r', sprintf('T_{t} = %.0f °C', Tt), 'LabelHorizontalAlignment','left');
xline(t_req, '--k', sprintf('t = %.1f s (%.2f min)', t_req, t_req/60), 'LabelOrientation','horizontal');
plot(t_req, Tt, 'ok', 'MarkerFaceColor','k');
xlabel('Tiempo (s)')
ylabel('Temperatura (°C)')

%% Ejemplo 6.4

% Conduccion transitoria en un solido semi-infinito expuesto a flujo constante
clear; clc;
% Datos
q_flux = 1250;               % [W/m2]
Ti = 20;                     % [°C]
t = 20*60;                   % [s]
    
% Estructura: nombre, conductividad k [W/m·K], difusividad alpha [m2/s]
materials = {
    'M 1',   23,   1.2e-5;
    'M 2', 265,    6.81e-5;
    'M 3',    38,     1.3e-5;
    'M 4',    574,    1.11e-4};
    
x = linspace(0,0.5,200);     % profundidad hasta 0.5 m
    
figure; hold on;
colors = {'-c','-k','-g','-m'};   % colores para la grafica
    
for i = 1:size(materials,1)
    name   = materials{i,1};
    k      = materials{i,2};
    alpha  = materials{i,3};
     
    % Ecuacion general (Caso 2, flujo constante en la superficie)
    T = Ti + ...
       (2*q_flux/k)*sqrt(alpha*t/pi).*exp(-x.^2./(4*alpha*t)) ...
       - (q_flux.*x./k).*erfc(x./(2*sqrt(alpha*t)));
    
   Tsup = T(1);
    
   fprintf('Ts en 20 min (%s): %.2f °C\n', ...
       name, Tsup);
      
        % Graficar
    plot(x,T,colors{i},'LineWidth',2,'DisplayName',name);
end
   
xlabel('Distancia x desde la superficie [m]');
ylabel('T [°C]');
title('Distribución de temperatura en sólido semi-infinito (t = 20 min)');
legend('show','Location','northeast');
grid on;

%% Ejemplo 6.5

function barra_cobre_explicito()
clc; clear; close all;
    
k = 400.0;
rho = 9000.0;
cp = 385.0;
alpha = k / (rho * cp);
L = 0.3;
Nx = 8;
dx = L / (Nx - 1);
x = linspace(0, L, Nx);
    
dt = 0.05;
t_final = 60.0;
nt = floor(t_final / dt);
    
r = alpha * dt / dx^2;
    
A = zeros(Nx, Nx);
for i = 2:(Nx-1)
   A(i,i) = 1 - 2*r;
   A(i,i-1) = r;
   A(i,i+1) = r;
end
    
A(1,1) = 1;
A(Nx,Nx) = 1;
    
T_left_init = 30.0;
T_right = 50.0;
T_left_new = 502.0;
    
T = T_left_init + (T_right - T_left_init) / L * x;
    
times_save = [0, 5, 30, 60];
num_saves = length(times_save);
saved = cell(1, num_saves);
saved{1} = T;
saved_flag = false(1, num_saves);
saved_flag(1) = true;
   
for n = 1:nt
    T = (A * T')'; 
      
    T(1) = T_left_new;
    T(Nx) = T_right;
        
    current_t = n * dt;
    for k_save = 1:num_saves
        ts = times_save(k_save);
            
            % Condición
        if abs(current_t - ts) < dt/2 && ~saved_flag(k_save) 
             saved{k_save} = T;
             saved_flag(k_save) = true;
        end
    end
end
    
if abs(nt * dt - t_final) < dt/2 && ~saved_flag(end)
    saved{end} = T;
end
    
% Tablas
Headers = cell(1, num_saves);
for i = 1:num_saves
   Headers{i} = sprintf('t=%.0fs', times_save(i));
end
    
Data_matrix = zeros(Nx, num_saves);
for i = 1:num_saves
    Data_matrix(:, i) = saved{i}';
end
    
disp(newline);
disp('Resultados de temperatura en nodos (columnas son tiempos):');
disp(Headers);
disp(round(Data_matrix, 2));
 
figure('Name', 'Evolución de T(x,t) en la barra de cobre (Explícito)', 'NumberTitle', 'off');
hold on;
Colores = lines(num_saves);
 
for i = 1:num_saves
   plot(x, saved{i}, 'Marker', 'o', 'Color', Colores(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('t=%.0fs', times_save(i)));
end
   
xlabel('x [m]');
ylabel('T [°C]');
title('Evolución de T(x,t) en la barra de cobre (Explícito)');
legend('Location', 'best');
grid on;
hold off;
    
end

%% Ejemplo 6.6.

function pared_implicta_transitorio()
    
clc; clear; close all;
    
%Datos
e = 0.65;
k = 18.0;
Cp = 0.12;
rho = 490.0;
sigma = 0.173e-8;
    
Ts0 = 1100.0;
Tinf = 90.0;
F_R = @(T_F) T_F + 459.67;
  
L = 8/12;
N = 21;
dx = L/(N-1);
   
dt = 0.01;
alpha = k/(rho*Cp);
Fo = alpha*dt/dx^2;
    
fprintf('Fo = %.4f\n', Fo);
  
hc = (Ts0 - Tinf)^0.25;
hr = sigma*e*((F_R(Ts0)^4 - F_R(Tinf)^4)/(Ts0 - Tinf));
    
Bi_c = hc*(dx/k);
Bi_r = hr*(dx/k);
   
A = zeros(N,N);
    
for i = 2:(N-1) 
    A(i,i-1) = -Fo;
    A(i,i)   = 1 + 2*Fo;
    A(i,i+1) = -Fo;
end
    
A(1,1) = 1 + 2*Fo;
A(1,2) = -2*Fo;
    
A(N,N-1) = -2*Fo;
A(N,N) = 1 + 2*Fo + 2*Fo*Bi_c + 2*Fo*Bi_r;
    
Ainv = inv(A);
T = ones(1, N) * Ts0; 
t_final = 4;
nt = floor(t_final/dt);
  
T_record = [];
time = [];
  
for n = 0:nt
   if mod(n, 200) == 0 
       T_record = [T_record; T];
       time = [time, n*dt];
    end
       
    b = T'; 
    b(N) = b(N) + (2*Fo*Bi_c + 2*Fo*Bi_r)*Tinf;
    T_nueva = Ainv * b;
    T = T_nueva';
      
end
    
if mod(nt, 200) ~= 0 
   T_record = [T_record; T];
   time = [time, nt*dt];
end
  
   
% Mapa de calor
figure('Name', 'Mapa de calor T(x,t)', 'NumberTitle', 'off', 'Position', [100 100 700 500]);
x_vals = linspace(0, L, N);
  
[X, Y] = meshgrid(x_vals, time);
    
contourf(X, Y, T_record, 50, 'LineStyle', 'none');
h = colorbar;
ylabel(h, 'Temperatura [°F]');
xlabel('Posición x [ft]');
ylabel('Tiempo [h]');
title('Mapa de calor T(x,t) (Contorno)');
colormap('jet'); 
axis tight;
    
end
