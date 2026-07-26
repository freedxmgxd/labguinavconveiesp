clear; clc; close all;

%% 1. Parâmetros do Problema
M = 1.0;            % Massa do módulo lunar (kg)
g = 2.0;            % Aceleração da gravidade lunar (m/s^2)
E = 3.0;            % Empuxo do retrofoguete (N)
theta_fix = 34;     % Ângulo limite (graus)

% Condições Iniciais: y0 = 70m, V0 = 4i - 7j
y0 = 70.0; V0x = 4.0; V0y = -7.0; x0 = 0.0;

%% 2. Acelerações (2ª Lei de Newton: m*ay = E*cos(theta) - m*g)
ax_neg  = (E * sind(-theta_fix)) / M;   % ax para theta = -34° (-1.6776 m/s^2)
ay_zero = (E * cosd(0)) / M - g;        % ay para theta = 0°  (1.0000 m/s^2)
ay_34   = (E * cosd(theta_fix)) / M - g;% ay para theta = 34° (0.4871 m/s^2)

%% 3. Janela de Pouso Suave (Torricelli: y = Vy^2 / (2*ay))
Vy_vec = linspace(-15, 0, 300);
y_parabola1 = (Vy_vec.^2) / (2 * ay_zero);  % Parábola 1: theta = 0°
y_parabola2 = (Vy_vec.^2) / (2 * ay_34);    % Parábola 2: theta = ±34°

figure('Name', 'Janela de Pouso Suave');
plot(Vy_vec, y_parabola1, 'r-', 'LineWidth', 2); hold on; grid on; box on;
plot(Vy_vec, y_parabola2, 'b--', 'LineWidth', 2);
plot(V0y, y0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
xlabel('V_y (m/s)'); ylabel('y (m)'); title('Janela de Pouso Suave');
legend('Parábola 1: \theta = 0°', 'Parábola 2: \theta = \pm 34°', 'C.I. (-7, 70)', 'Location', 'northwest');

%% 4. Estratégias de Pouso Suave
t_freio_x = (0 - V0x) / ax_neg;         % Tempo a -34° para zerar Vx (2.38 s)
dx_brake  = V0x * t_freio_x + 0.5 * ax_neg * t_freio_x^2;

V2_base = V0y + ay_34 * t_freio_x;
y_base  = y0 + V0y * t_freio_x + 0.5 * ay_34 * t_freio_x^2;

% --- Estratégia 1: Disparo Único -> Parábola 1 (\theta = 0°) ---
A1 = (g/2) * (g/ay_zero + 1);
B1 = -(g * V2_base)/ay_zero - (V0y - g * t_freio_x);
C1 = (V2_base^2)/(2 * ay_zero) - y_base;

t_livre1 = (-B1 + sqrt(B1^2 - 4*A1*C1)) / (2*A1);
Vy_ch1   = V0y - g*t_livre1 + ay_34*t_freio_x;
y_ch1    = Vy_ch1^2 / (2*ay_zero);
x_ch1    = V0x * t_livre1 + dx_brake;
t_freio1_y = (0 - Vy_ch1) / ay_zero;
t1_tot   = t_livre1 + t_freio_x + t_freio1_y;
t1_queima= t_freio_x + t_freio1_y;

% --- Estratégia 2: Disparo Único -> Parábola 2 (\theta = ±34°) ---
A2 = (g/2) * (g/ay_34 + 1);
B2 = -(g * V2_base)/ay_34 - (V0y - g * t_freio_x);
C2 = (V2_base^2)/(2 * ay_34) - y_base;

t_livre2 = (-B2 + sqrt(B2^2 - 4*A2*C2)) / (2*A2);
Vy_ch2   = V0y - g*t_livre2 + ay_34*t_freio_x;
y_ch2    = Vy_ch2^2 / (2*ay_34);
x_ch2    = V0x * t_livre2 + dx_brake;
t_freio2_y = (0 - Vy_ch2) / ay_34;
t2_tot   = t_livre2 + t_freio_x + t_freio2_y;
t2_queima= t_freio_x + t_freio2_y;

% --- Estratégia 3: Disparo Misto (\theta = ±34° -> 0°) ---
t_livre3 = (t_livre1 + t_livre2) / 2.0;
Vy_ch3   = V0y - g*t_livre3 + ay_34*t_freio_x;
y_ch3    = (y0 + V0y*t_livre3 - 0.5*g*t_livre3^2) + (V0y - g*t_livre3)*t_freio_x + 0.5*ay_34*t_freio_x^2;
x_ch3    = V0x * t_livre3 + dx_brake;

A3_quad = 0.5 * ay_34 * (1 - ay_34 / ay_zero);
B3_quad = Vy_ch3 * (1 - ay_34 / ay_zero);
C3_quad = y_ch3 - (Vy_ch3^2) / (2 * ay_zero);

t3_A = (-B3_quad - sqrt(B3_quad^2 - 4*A3_quad*C3_quad)) / (2*A3_quad);
t3_B = (0 - Vy_ch3 - ay_34 * t3_A) / ay_zero;
t3_queima = t_freio_x + t3_A + t3_B;
t3_tot   = t_livre3 + t3_queima;

%% 5. Impressão dos Resultados
fprintf('=================================================================\n');
fprintf('            ANÁLISE DAS ESTRATÉGIAS DE POUSO SUAVE               \n');
fprintf('=================================================================\n\n');

fprintf('Estratégia 1 (Disparo Único -> Parábola 1: \\theta = 0°):\n');
fprintf('  1) t = 0.00s a %.2fs : Foguete DESLIGADO (Queda livre inicial)\n', t_livre1);
fprintf('  2) t = %.2fs a %.2fs : \\theta = -34° (LIGA FOGUETE: freia Vx)\n', t_livre1, t_livre1 + t_freio_x);
fprintf('  3) t = %.2fs a %.2fs : \\theta = 0°   (Motor LIGADO sobre Parábola 1)\n', t_livre1 + t_freio_x, t1_tot);
fprintf('  -> Tempo Total: %.2f s | Tempo de Queima: %.2f s | Posição Final xf: %.2f m\n\n', t1_tot, t1_queima, x_ch1);

fprintf('Estratégia 2 (Disparo Único -> Parábola 2: \\theta = \\pm34°):\n');
fprintf('  1) t = 0.00s a %.2fs : Foguete DESLIGADO (Queda livre inicial)\n', t_livre2);
fprintf('  2) t = %.2fs a %.2fs : \\theta = -34° (LIGA FOGUETE: freia Vx)\n', t_livre2, t_livre2 + t_freio_x);
fprintf('  3) t = %.2fs a %.2fs : \\theta = \\pm34° (Motor LIGADO sobre Parábola 2)\n', t_livre2 + t_freio_x, t2_tot);
fprintf('  -> Tempo Total: %.2f s | Tempo de Queima: %.2f s | Posição Final xf: %.2f m\n\n', t2_tot, t2_queima, x_ch2);

fprintf('Estratégia 3 (Disparo Misto -> Transição \\pm34° -> 0°):\n');
fprintf('  1) t = 0.00s a %.2fs : Foguete DESLIGADO (Queda livre inicial)\n', t_livre3);
fprintf('  2) t = %.2fs a %.2fs : \\theta = -34° (LIGA FOGUETE: freia Vx)\n', t_livre3, t_livre3 + t_freio_x);
fprintf('  3) t = %.2fs a %.2fs : \\theta = \\pm34° (3.70s) e \\theta = 0° (5.92s) até pouso\n', t_livre3 + t_freio_x, t3_tot);
fprintf('  -> Tempo Total: %.2f s | Tempo de Queima: %.2f s | Posição Final xf: %.2f m\n\n', t3_tot, t3_queima, x_ch3);

fprintf('=================================================================\n');
fprintf('CONCLUSÃO (Questão 6):\n');
fprintf('A ESTRATÉGIA 1 É A MELHOR E MINIMIZA O GASTO DE COMBUSTÍVEL!\n');
fprintf('  - Tempo de queima Est. 1: %.2f s (Vencedora)\n', t1_queima);
fprintf('  - Tempo de queima Est. 3: %.2f s\n', t3_queima);
fprintf('  - Tempo de queima Est. 2: %.2f s\n', t2_queima);
fprintf('A Estratégia 1 economiza %.2f s de queima em relação à Estratégia 2\n', t2_queima - t1_queima);
fprintf('e %.2f s em relação à Estratégia 3, pois o empuxo a \\theta = 0° maximiza\n', t3_queima - t1_queima);
fprintf('a aceleração vertical (ay = 1.0 m/s^2).\n');
fprintf('=================================================================\n\n');

%% 6. Simulação e Gráficos
dt = 0.005;

% Simulação Estratégia 1
t_sim1 = 0:dt:t1_tot; N1 = length(t_sim1);
x1 = zeros(1,N1); y1_sim = zeros(1,N1); Vx1 = zeros(1,N1); Vy1_sim = zeros(1,N1);
x1(1) = x0; y1_sim(1) = y0; Vx1(1) = V0x; Vy1_sim(1) = V0y;
for k = 1:N1-1
    t = t_sim1(k);
    if t <= t_livre1, theta = 0; E_c = 0;
    elseif t <= (t_livre1 + t_freio_x), theta = -theta_fix; E_c = E;
    else, theta = 0; E_c = E; end
    ax = (E_c * sind(theta))/M; ay = (E_c * cosd(theta))/M - g;
    Vx1(k+1) = Vx1(k) + ax*dt; Vy1_sim(k+1) = Vy1_sim(k) + ay*dt;
    x1(k+1)  = x1(k) + Vx1(k)*dt; y1_sim(k+1) = max(0, y1_sim(k) + Vy1_sim(k)*dt);
end

% Simulação Estratégia 2
t_sim2 = 0:dt:t2_tot; N2 = length(t_sim2);
x2 = zeros(1,N2); y2_sim = zeros(1,N2); Vx2 = zeros(1,N2); Vy2_sim = zeros(1,N2);
x2(1) = x0; y2_sim(1) = y0; Vx2(1) = V0x; Vy2_sim(1) = V0y;
for k = 1:N2-1
    t = t_sim2(k);
    if t <= t_livre2, theta = 0; E_c = 0;
    elseif t <= (t_livre2 + t_freio_x), theta = -theta_fix; E_c = E;
    else, theta = theta_fix; E_c = E; end
    if t > (t_livre2 + t_freio_x), ax = 0; else, ax = (E_c * sind(theta))/M; end
    ay = (E_c * cosd(theta))/M - g;
    Vx2(k+1) = Vx2(k) + ax*dt; Vy2_sim(k+1) = Vy2_sim(k) + ay*dt;
    x2(k+1)  = x2(k) + Vx2(k)*dt; y2_sim(k+1) = max(0, y2_sim(k) + Vy2_sim(k)*dt);
end

% Simulação Estratégia 3
t_sim3 = 0:dt:t3_tot; N3 = length(t_sim3);
x3 = zeros(1,N3); y3_sim = zeros(1,N3); Vx3 = zeros(1,N3); Vy3_sim = zeros(1,N3);
x3(1) = x0; y3_sim(1) = y0; Vx3(1) = V0x; Vy3_sim(1) = V0y;
for k = 1:N3-1
    t = t_sim3(k);
    if t <= t_livre3, theta = 0; E_c = 0;
    elseif t <= (t_livre3 + t_freio_x), theta = -theta_fix; E_c = E;
    elseif t <= (t_livre3 + t_freio_x + t3_A), theta = theta_fix; E_c = E;
    else, theta = 0; E_c = E; end
    if t > (t_livre3 + t_freio_x), ax = 0; else, ax = (E_c * sind(theta))/M; end
    ay = (E_c * cosd(theta))/M - g;
    Vx3(k+1) = Vx3(k) + ax*dt; Vy3_sim(k+1) = Vy3_sim(k) + ay*dt;
    x3(k+1)  = x3(k) + Vx3(k)*dt; y3_sim(k+1) = max(0, y3_sim(k) + Vy3_sim(k)*dt);
end

% Figura 2: Estratégia 1
figure('Name', 'Estratégia 1', 'Color', 'w');
subplot(1,2,1); plot(x1, y1_sim, 'b-', 'LineWidth', 2.5); hold on; grid on; box on;
plot(x1(1), y1_sim(1), 'ko', 'MarkerFaceColor', 'k');
plot(x1(end), y1_sim(end), 'r*', 'MarkerSize', 10);
xlabel('x (m)'); ylabel('y (m)'); title('Estratégia 1: Trajetória Espacial');
legend('Trajetória', 'C.I. (0, 70)', 'Pouso (yf=0)');

subplot(1,2,2); plot(Vy_vec, y_parabola1, 'r--', 'LineWidth', 1.5); hold on; grid on; box on;
plot(Vy_vec, y_parabola2, 'b--', 'LineWidth', 1.5);
plot(Vy1_sim, y1_sim, 'b-', 'LineWidth', 2.5); plot(V0y, y0, 'ko', 'MarkerFaceColor', 'k');
xlabel('V_y (m/s)'); ylabel('y (m)'); title('Estratégia 1: Plano de Fase');
legend('Parábola 1', 'Parábola 2', 'Est 1', 'C.I.', 'Location', 'northwest');

% Figura 3: Estratégia 2
figure('Name', 'Estratégia 2', 'Color', 'w');
subplot(1,2,1); plot(x2, y2_sim, 'g-', 'LineWidth', 2.5); hold on; grid on; box on;
plot(x2(1), y2_sim(1), 'ko', 'MarkerFaceColor', 'k');
plot(x2(end), y2_sim(end), 'r*', 'MarkerSize', 10);
xlabel('x (m)'); ylabel('y (m)'); title('Estratégia 2: Trajetória Espacial');
legend('Trajetória', 'C.I. (0, 70)', 'Pouso (yf=0)');

subplot(1,2,2); plot(Vy_vec, y_parabola1, 'r--', 'LineWidth', 1.5); hold on; grid on; box on;
plot(Vy_vec, y_parabola2, 'b--', 'LineWidth', 1.5);
plot(Vy2_sim, y2_sim, 'g-', 'LineWidth', 2.5); plot(V0y, y0, 'ko', 'MarkerFaceColor', 'k');
xlabel('V_y (m/s)'); ylabel('y (m)'); title('Estratégia 2: Plano de Fase');
legend('Parábola 1', 'Parábola 2', 'Est 2', 'C.I.', 'Location', 'northwest');

% Figura 4: Estratégia 3
figure('Name', 'Estratégia 3', 'Color', 'w');
subplot(1,2,1); plot(x3, y3_sim, 'm-', 'LineWidth', 2.5); hold on; grid on; box on;
plot(x3(1), y3_sim(1), 'ko', 'MarkerFaceColor', 'k');
plot(x3(end), y3_sim(end), 'r*', 'MarkerSize', 10);
xlabel('x (m)'); ylabel('y (m)'); title('Estratégia 3: Trajetória Espacial');
legend('Trajetória', 'C.I. (0, 70)', 'Pouso (yf=0)');

subplot(1,2,2); plot(Vy_vec, y_parabola1, 'r--', 'LineWidth', 1.5); hold on; grid on; box on;
plot(Vy_vec, y_parabola2, 'b--', 'LineWidth', 1.5);
plot(Vy3_sim, y3_sim, 'm-', 'LineWidth', 2.5); plot(V0y, y0, 'ko', 'MarkerFaceColor', 'k');
xlabel('V_y (m/s)'); ylabel('y (m)'); title('Estratégia 3: Plano de Fase');
legend('Parábola 1', 'Parábola 2', 'Est 3', 'C.I.', 'Location', 'northwest');

saveas(figure(1), 'janela_pouso.png');
saveas(figure(2), 'estrategia1.png');
saveas(figure(3), 'estrategia2.png');
saveas(figure(4), 'estrategia3.png');
