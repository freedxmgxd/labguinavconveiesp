clear; clc; close all;

% =========================================================================
% 1. DADOS DE ENTRADA E CONSTANTES
% =========================================================================

% TODO: revisar - dados abaixo eram de rascunho/teste anterior (verificar origem)
% r = [8228 389 6888];      % Posição inicial (km)
% v = [-0.7 6.6 -0.6];      % Velocidade inicial (km/s)

% Dados do enunciado (problema2-gnc.profGil.pdf)
r = [-375.12842 -5011.13118 4469.62763];  % Posição inicial (km)
v = [7.471975   -1.529935   -1.077989];   % Velocidade inicial (km/s)

mi          = 3.986e5;          % Parâmetro gravitacional da Terra (km^3/s^2)
Raio_T      = 6378.137;         % Raio equatorial da Terra (km)
omega_terra = 7.292115e-5;      % Velocidade angular da Terra (rad/s)

norm_r = norm(r);
norm_v = norm(v);

% =========================================================================
% 2. ELEMENTOS ORBITAIS CLÁSSICOS (EOC / COE)
% =========================================================================
epsilon = (norm_v^2 / 2) - (mi / norm_r);   % Energia específica (km^2/s^2)
a       = -mi / (2 * epsilon);               % Semieixo maior (km)

exc       = (1/mi) * ((norm_v^2 - mi/norm_r) * r - dot(r, v) * v);  % Vetor excentricidade
norma_exc = norm(exc);                       % Excentricidade (adimensional)

h     = cross(r, v);                         % Momento angular específico (km^2/s)
h_mag = norm(h);
p     = h_mag^2 / mi;                        % Semi-latus rectum (km)

K     = [0 0 1];
inc   = acosd(dot(K, h) / h_mag);           % Inclinação (graus)

n     = cross(K, h);                         % Vetor nó ascendente
i_vec = [1 0 0];

% RAAN (Ascensão Reta do Nó Ascendente)
ARNA = acosd(dot(i_vec, n) / norm(n));
if n(2) < 0, ARNA = 360.0 - ARNA; end

% Argumento do Perigeu (omega)
OMEGA = acosd(dot(n, exc) / (norm(n) * norma_exc));
if exc(3) < 0, OMEGA = 360.0 - OMEGA; end

% Anomalia Verdadeira inicial
NI = acosd(dot(exc, r) / (norma_exc * norm_r));
if dot(r, v) < 0, NI = 360.0 - NI; end

% Grandezas derivadas
rp  = p / (1 + norma_exc);             % Raio do perigeu (km)
ra  = p / (1 - norma_exc);             % Raio do apogeu (km)
n_medio = sqrt(mi / a^3);              % Movimento médio (rad/s)
T   = 2 * pi / n_medio;               % Período orbital (s)

% =========================================================================
% 3. IMPRESSÃO DOS EOC NO COMMAND WINDOW (Item 1 do enunciado)
% =========================================================================
fprintf('\n==========================================================\n');
fprintf('       ELEMENTOS ORBITAIS CLÁSSICOS (EOC)\n');
fprintf('==========================================================\n');
fprintf('  Semieixo maior           a   = %12.4f  km\n',   a);
fprintf('  Excentricidade           e   = %12.6f\n',       norma_exc);
fprintf('  Inclinação               i   = %12.4f  graus\n',inc);
fprintf('  RAAN (Omega maiúsculo)   Ω   = %12.4f  graus\n',ARNA);
fprintf('  Arg. do Perigeu          ω   = %12.4f  graus\n',OMEGA);
fprintf('  Anomalia Verdadeira ini. ν₀  = %12.4f  graus\n',NI);
fprintf('----------------------------------------------------------\n');
fprintf('  Semi-latus rectum        p   = %12.4f  km\n',   p);
fprintf('  Raio do perigeu          rp  = %12.4f  km\n',   rp);
fprintf('  Raio do apogeu           ra  = %12.4f  km\n',   ra);
fprintf('  Alt. perigeu (aprox.)        = %12.4f  km\n',   rp - Raio_T);
fprintf('  Alt. apogeu  (aprox.)        = %12.4f  km\n',   ra - Raio_T);
fprintf('  Período orbital          T   = %12.2f  s\n',    T);
fprintf('  Momento angular          |h| = %12.4f  km²/s\n',h_mag);
fprintf('==========================================================\n\n');

% =========================================================================
% 4. GERAÇÃO DA ÓRBITA NO PLANO PERIFOCAL
% =========================================================================
ni_vetor = linspace(0, 2*pi, 1000);
r_polar  = p ./ (1 + norma_exc * cos(ni_vetor));

% Coordenadas no plano orbital (X_perifocal, Y_perifocal)
Ppolar = [r_polar .* cos(ni_vetor); r_polar .* sin(ni_vetor); zeros(1, length(ni_vetor))];

% Pontos especiais: perigeu (ni=0) e apogeu (ni=pi)
r_perigeu = p / (1 + norma_exc);
r_apogeu  = p / (1 - norma_exc);
P_perigeu = [r_perigeu; 0; 0];
P_apogeu  = [-r_apogeu; 0; 0];

% =========================================================================
% 5. ROTAÇÃO DO PLANO PERIFOCAL PARA O SGI (Inercial)
% =========================================================================
% Convenção conforme PB2-plano.orbital-SGI.pdf (Prof. Gil):
%   Rotações no sentido HORÁRIO => ângulos NEGATIVOS.
%   1ª: -omega em torno de z  (rotz1)
%   2ª: -i     em torno de x  (rotx2)
%   3ª: -arna  em torno de z  (rotz3)
%   Pxyz = rotz3 * rotx2 * rotz1 * Ppolar
rotz1 = [cosd(-OMEGA)  sind(-OMEGA) 0; -sind(-OMEGA) cosd(-OMEGA) 0; 0 0 1];
rotx2 = [1 0 0; 0 cosd(-inc)  sind(-inc); 0 -sind(-inc) cosd(-inc)];
rotz3 = [cosd(-ARNA)   sind(-ARNA)  0; -sind(-ARNA)  cosd(-ARNA)  0; 0 0 1];
Q     = rotz3 * rotx2 * rotz1;

% Órbita completa em SGI
Pxyz = Q * Ppolar;
rx   = Pxyz(1, :);
ry   = Pxyz(2, :);
rz   = Pxyz(3, :);

% Perigeu e apogeu em SGI
sgi_perigeu = Q * P_perigeu;
sgi_apogeu  = Q * P_apogeu;

% ----- Verificação numérica do fechamento da órbita -----
% O ponto r (dado do enunciado) deve coincidir com o ponto da órbita em NI.
% Erro < 1 km é aceitável dado o passo discreto de ni_vetor.
idx_ni = max(1, round(NI / 360 * (length(ni_vetor)-1)) + 1);
erro_fechamento = norm([rx(idx_ni) ry(idx_ni) rz(idx_ni)] - r);
fprintf('\n[VERIFICAÇÃO] Erro de fechamento da órbita em NI=%.4f°: %.4f km\n', NI, erro_fechamento);
if erro_fechamento > 10
    warning('Erro de fechamento alto (%.2f km). Verifique as rotações!', erro_fechamento);
end

% =========================================================================
% 6. CÁLCULO DE LATITUDE E LONGITUDE (C/ ROTAÇÃO DA TERRA)
% =========================================================================
lat = atan2d(rz, sqrt(rx.^2 + ry.^2));

% Tempo de voo desde o perigeu para cada ponto (Equação de Kepler)
E_vetor = 2 * atan(sqrt((1 - norma_exc)/(1 + norma_exc)) * tan(ni_vetor/2));
M_vetor = E_vetor - norma_exc * sin(E_vetor);
t       = M_vetor / n_medio;
t(t < 0) = t(t < 0) + 2*pi/n_medio;   % Ajuste de tempo negativo

% TODO: revisar - longitude estava em [0, 360]; ajustado para [-180, 180] para
%                 compatibilidade com mapas padrão (mapa-mundi)
% lon = atan2d(ry, rx) - rad2deg(omega_terra * t);
% lon = mod(lon, 360);   % TODO: revisar - mod para [0,360] causava artefato visual

lon = atan2d(ry, rx) - rad2deg(omega_terra * t);
lon = mod(lon + 180, 360) - 180;   % Intervalo [-180°, 180°]

% Tratar descontinuidades de longitude para plotagem contínua
% Inserir NaN onde o salto é maior que 180° (wrap-around do mapa)
dlon       = abs(diff(lon));
idx_breaks = find(dlon > 180);
lon_plot   = lon;
lat_plot   = lat;
lon_plot(idx_breaks) = NaN;
lat_plot(idx_breaks) = NaN;

% =========================================================================
% 7. GRÁFICOS
% =========================================================================

% --- Gráfico 1: Plano 2D Perifocal ---
figure(1);
plot(Ppolar(1,:)/1e3, Ppolar(2,:)/1e3, 'b-', 'LineWidth', 1.8); hold on; grid on; axis equal;
plot(0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');           % Terra (foco)
plot(P_perigeu(1)/1e3, P_perigeu(2)/1e3, 'g^', 'MarkerSize', 9, 'MarkerFaceColor', 'g'); % Perigeu
plot(P_apogeu(1)/1e3,  P_apogeu(2)/1e3,  'ms',  'MarkerSize', 9, 'MarkerFaceColor', 'm'); % Apogeu
legend('Órbita', 'Terra (foco)', 'Perigeu', 'Apogeu', 'Location', 'best');
title('Plano Perifocal 2D');
xlabel('X_p (10³ km)'); ylabel('Y_p (10³ km)');
hold off;

% --- Gráfico 2: Órbita em 3D (SGI) ---
figure(2);
plot3(rx, ry, rz, 'b-', 'LineWidth', 1.8); hold on; grid on; axis equal; view(3);
plot3(r(1), r(2), r(3), 'mo', 'MarkerSize', 9, 'MarkerFaceColor', 'm');          % Posição inicial
plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');                  % Terra
plot3(sgi_perigeu(1), sgi_perigeu(2), sgi_perigeu(3), 'g^', ...
      'MarkerSize', 9, 'MarkerFaceColor', 'g');                                   % Perigeu
plot3(sgi_apogeu(1),  sgi_apogeu(2),  sgi_apogeu(3),  'ms', ...
      'MarkerSize', 9, 'MarkerFaceColor', 'm');                                    % Apogeu
legend('Órbita (SGI)', 'Pos. inicial', 'Terra', 'Perigeu', 'Apogeu', 'Location', 'best');
title('Sistema Geocêntrico Inercial 3D (Trajetória)');
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
hold off;

% --- Gráfico 3: Rascunho do Sistema Cartesiano Terrestre (Fig. B.2 do Apêndice) ---
figure(3);
hold on; grid on; axis equal; view(3);

% Esfera simplificada para representar a Terra
[X_esf, Y_esf, Z_esf] = sphere(15);
mesh(X_esf*Raio_T, Y_esf*Raio_T, Z_esf*Raio_T, 'EdgeColor', [0.8 0.8 0.8], 'FaceColor', 'none');

% Eixos coordenados fixos
max_ax = 9000;
plot3([0 max_ax], [0 0], [0 0], 'k', 'LineWidth', 1.2);
plot3([0 0], [0 max_ax], [0 0], 'k', 'LineWidth', 1.2);
plot3([0 0], [0 0], [0 max_ax], 'k', 'LineWidth', 1.2);
text(max_ax+300, 0, 0, 'X_T');
text(0, max_ax+300, 0, 'Y_T');
text(0, 0, max_ax+300, 'Z_T');

% Ponto inicial P_int e vetor posição
px = r(1); py = r(2); pz = r(3);
plot3([0 px], [0 py], [0 pz], 'k-', 'LineWidth', 1.5);
plot3(px, py, pz, 'rx', 'MarkerSize', 8, 'LineWidth', 2);
text(px+300, py, pz+300, 'P_{int}');

% Linhas de projeção geométrica
plot3([0 px], [0 py], [0 0], 'g-',  'LineWidth', 1.2);   % Projeção no plano XY
plot3([px px], [py py], [0 pz], 'g--', 'LineWidth', 1);  % Componente Z
plot3([px px], [0 py], [0 0], 'g:',  'LineWidth', 1);    % Linha paralela a Y
plot3([0 px], [py py], [0 0], 'g:',  'LineWidth', 1);    % Linha paralela a X

% Arcos indicativos de Latitude e Longitude
r_arco         = 2500;
lambda_E_rad   = atan2(py, px);
if lambda_E_rad < 0, lambda_E_rad = 2*pi + lambda_E_rad; end

theta_arc = linspace(0, lambda_E_rad, 20);
plot3(r_arco*cos(theta_arc), r_arco*sin(theta_arc), zeros(1,20), 'b-', 'LineWidth', 1.5);
text(r_arco*cos(lambda_E_rad/2)+200, r_arco*sin(lambda_E_rad/2)+200, 0, '\lambda_E', ...
     'Color', 'b', 'FontWeight', 'bold');

r_xy        = sqrt(px^2 + py^2);
phi_linha_rad = atan2(pz, r_xy);
phi_arc     = linspace(0, phi_linha_rad, 20);
plot3(r_arco*cos(phi_arc)*px/r_xy, r_arco*cos(phi_arc)*py/r_xy, r_arco*sin(phi_arc), ...
      'r-', 'LineWidth', 1.5);
text(r_arco*px/r_xy, r_arco*py/r_xy, r_arco*sin(phi_linha_rad/2)+300, '\Phi''', ...
     'Color', 'r', 'FontWeight', 'bold');

title('Fig B.2 - Rascunho (Latitude e Longitude Geocêntricas)');
hold off;

% --- Gráfico 4: Rastreio no Solo (Ground Track) ---
figure(4);

% Mapa de fundo: linhas de costa (disponível no MATLAB base)
if exist('coastlines', 'file') || exist('coast', 'file')
    try
        load coastlines;           % Variáveis: coastlat, coastlon
        plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);
    catch
        % fallback: sem mapa
    end
else
    % Alternativa sem Mapping Toolbox: grade simples de referência
    [lon_g, lat_g] = meshgrid(-180:30:180, -90:30:90);
    plot(lon_g', lat_g', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5); hold on;
    plot(lon_g,  lat_g,  'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
end
hold on; grid on;

% Ground track (linha contínua, com NaN nas descontinuidades)
% TODO: revisar - antes era plot com pontos '.'; agora usa linha com NaN nos saltos
% plot(lon, lat, 'b.', 'MarkerSize', 6);
plot(lon_plot, lat_plot, 'b-', 'LineWidth', 1.4);

% Ponto inicial
plot(lon(1), lat(1), 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
legend('Linhas de costa', 'Ground Track', 'Pos. inicial', 'Location', 'best');

title('Rastreio no Solo Real (Ground Track)');
xlabel('Longitude (graus)'); ylabel('Latitude (graus)');

% TODO: revisar - eixo estava em [0, 360]; ajustado para [-180, 180]
% axis([0 360 -90 90]); xticks(0:45:360);
axis([-180 180 -90 90]);
xticks(-180:45:180);
yticks(-90:30:90);
hold off;