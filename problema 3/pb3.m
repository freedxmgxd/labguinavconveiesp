function pb3()
    clear; clc; close all;

    % =========================================================================
    % 1. CONSTANTES E PARÂMETROS DE ENTRADA
    % =========================================================================
    fprintf('==========================================================\n');
    fprintf('       LABORATÓRIO DE GN&C - PROBLEMA 3\n');
    fprintf('==========================================================\n\n');

    % Coordenadas Geográficas da Estação (Xichang Satellite Launch Center)
    lat_deg = 28; lat_min = 14; lat_sec = 45.66; % Latitude Geodésica N
    lon_deg = 102; lon_min = 01; lon_sec = 35.6; % Longitude Leste E
    alt_m = 1800;                                % Altitude (metros)

    % Instante da Medição: 20/06/2026 às 21h43m35.0s UT
    ano = 2026; mes = 6; dia = 20;
    hora = 21; minuto = 43; segundo = 35.0;

    % Conversão das Coordenadas para Graus Decimais
    phi = lat_deg + lat_min/60 + lat_sec/3600;       % Latitude geodésica (graus)
    lambda_E = lon_deg + lon_min/60 + lon_sec/3600;  % Longitude Leste (graus)
    H = alt_m / 1000;                                % Altitude em km

    % Parâmetros do Elipsoide de Referência WGS-84
    ae = 6378.137;                  % Semieixo maior (km)
    f = 1 / 298.257223563;          % Achatamento
    e2 = 2*f - f^2;                 % Primeira excentricidade linear ao quadrado

    % =========================================================================
    % 2. CÁLCULO DE DATA JULIANA (JD) E DATA JULIANA MODIFICADA (DJM)
    % =========================================================================
    [JD_meiodia, DJM0, DJM_frac] = DiaJulianoLocal(ano, mes, dia, hora, minuto, segundo);
    JD_frac = DJM_frac + 2400000.5;

    % =========================================================================
    % 3. CÁLCULO DO TEMPO SIDERAL DE GREENWICH (TSG)
    % =========================================================================
    tsg_graus = TSGLocal(ano, mes, dia, hora, minuto, segundo);

    % =========================================================================
    % 4. SISTEMA DE COORDENADAS CARTESIANAS TERRESTRES (SCT)
    % =========================================================================
    % Cálculo de N (Grande Normal / Raio de curvatura no primeiro vertical)
    N = ae / sqrt(1 - e2 * sind(phi)^2);

    % Coordenadas Cartesianas no SCT (aT)
    xT = (N + H) * cosd(phi) * cosd(lambda_E);
    yT = (N + H) * cosd(phi) * sind(lambda_E);
    zT = (N * (1 - e2) + H) * sind(phi);
    aT = [xT; yT; zT];

    % =========================================================================
    % 5. ROTAÇÃO DO SCT PARA O REFERENCIAL GEOCÊNTRICO INERCIAL (SGI)
    % =========================================================================
    % Rotação anti-horária de +tsg em torno do eixo Z:
    % aI_pos = Rz(tsg) * aT
    Rz = [cosd(tsg_graus), -sind(tsg_graus), 0; ...
          sind(tsg_graus),  cosd(tsg_graus), 0; ...
          0,                0,               1];

    aI_pos = Rz * aT;
    aI = aI_pos / norm(aI_pos); % Versor normalizado

    % =========================================================================
    % 6. DETERMINAÇÃO DAS COORDENADAS DE APONTAMENTO CELESTE (AR e DEC)
    % =========================================================================
    ra_deg = mod(atan2d(aI(2), aI(1)), 360);
    dec_deg = asind(aI(3));

    % Formatação de AR (h m s)
    ra_hours = ra_deg / 15;
    ar_h = floor(ra_hours);
    ar_m = floor((ra_hours - ar_h) * 60);
    ar_s = ((ra_hours - ar_h) * 60 - ar_m) * 60;

    % Formatação de DEC (d m s)
    dec_d = floor(dec_deg);
    dec_m = floor((dec_deg - dec_d) * 60);
    dec_s = ((dec_deg - dec_d) * 60 - dec_m) * 60;

    % =========================================================================
    % IMPRESSÃO DOS RESULTADOS NO CONSOLE
    % =========================================================================
    fprintf('----------------------------------------------------------\n');
    fprintf('1. Dados Geográficos de Entrada:\n');
    fprintf('   Latitude Geodésica (phi)  : %02d° %02d'' %05.2f\" N  (= %12.8f°)\n', lat_deg, lat_min, lat_sec, phi);
    fprintf('   Longitude Leste (lambda_E): %02d° %02d'' %05.2f\" E  (= %12.8f°)\n', lon_deg, lon_min, lon_sec, lambda_E);
    fprintf('   Altitude (H)              : %6.1f m         (= %12.6f km)\n', alt_m, H);
    fprintf('----------------------------------------------------------\n');
    fprintf('2. Tempo e Data Juliana (UT: %02d/%02d/%04d %02d:%02d:%04.1f):\n', dia, mes, ano, hora, minuto, segundo);
    fprintf('   Data Juliana (meio-dia)   : %d\n', JD_meiodia);
    fprintf('   Data Juliana com fração   : %18.10f\n', JD_frac);
    fprintf('   Data Juliana Modificada   : %18.10f (MJD)\n', DJM_frac);
    fprintf('   Tempo Sideral de Greenwich: %12.8f° (TSG)\n', tsg_graus);
    fprintf('----------------------------------------------------------\n');
    fprintf('3. Vetor no Sistema de Coordenadas Terrestres (SCT):\n');
    fprintf('   aT = [%14.6f, %14.6f, %14.6f] km\n', aT(1), aT(2), aT(3));
    fprintf('   Magnitude |aT|: %10.4f km\n', norm(aT));
    fprintf('----------------------------------------------------------\n');
    fprintf('4. Vetor no Sistema Geocêntrico Inercial (SGI):\n');
    fprintf('   aI_pos = [%14.6f, %14.6f, %14.6f] km\n', aI_pos(1), aI_pos(2), aI_pos(3));
    fprintf('   Versor aI = [%12.8f, %12.8f, %12.8f]\n', aI(1), aI(2), aI(3));
    fprintf('   Verificação (|aI|): %1.8f\n', norm(aI));
    fprintf('----------------------------------------------------------\n');
    fprintf('5. Coordenadas Celestiais Zenitais do Local (Apontamento):\n');
    fprintf('   Ascensão Reta (AR)  : %02dh %02dm %05.2fs  (= %12.8f°)\n', ar_h, ar_m, ar_s, ra_deg);
    fprintf('   Declinação (DEC)    : %+02d° %02d'' %05.2f\"   (= %12.8f°)\n', dec_d, dec_m, dec_s, dec_deg);
    fprintf('   Constelação no Zênite: PEGASUS\n');
    fprintf('==========================================================\n\n');

    % =========================================================================
    % 7. VISUALIZAÇÃO GRÁFICA 3D
    % =========================================================================
    figure(1);
    hold on; grid on; axis equal; view(135, 25);

    % Plotagem da Esfera Terrestre (Raio Médio)
    [X_esf, Y_esf, Z_esf] = sphere(30);
    surf(X_esf * ae, Y_esf * ae, Z_esf * ae, 'FaceColor', [0.8, 0.9, 1.0], ...
         'EdgeColor', [0.6, 0.7, 0.8], 'FaceAlpha', 0.25);

    % Plotagem do Equador e do Meridiano de Greenwich
    t_circle = linspace(0, 2*pi, 100);
    plot3(ae*cos(t_circle), ae*sin(t_circle), zeros(1, 100), 'k--', 'LineWidth', 0.8); % Equador
    plot3(ae*cos(t_circle), zeros(1, 100), ae*sin(t_circle), 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':', 'LineWidth', 0.8); % Greenwich

    % Plotagem do Mapa Mundi (Linhas de Costa)
    if exist('coastlines.mat', 'file')
        load('coastlines.mat'); % Carrega 'tmp' contendo [longitude, latitude]
        lon_c = tmp(:, 1);
        lat_c = tmp(:, 2);
        N_c = ae ./ sqrt(1 - e2 * sind(lat_c).^2);
        x_c = N_c .* cosd(lat_c) .* cosd(lon_c);
        y_c = N_c .* cosd(lat_c) .* sind(lon_c);
        z_c = N_c .* (1 - e2) .* sind(lat_c);
        plot3(x_c, y_c, z_c, 'Color', [0.3, 0.4, 0.6], 'LineWidth', 0.8);
    end

    % Comprimento dos Eixos Coordenados para Visualização
    L_eixos = 9000;

    % Eixos SCT (Fixos na Terra - Cor Preta)
    plot3([0 L_eixos], [0 0], [0 0], 'k-', 'LineWidth', 1.5);
    plot3([0 0], [0 L_eixos], [0 0], 'k-', 'LineWidth', 1.5);
    plot3([0 0], [0 0], [0 L_eixos], 'k-', 'LineWidth', 1.5);
    text(L_eixos + 300, 0, 0, 'X_T (Greenwich)', 'FontWeight', 'bold');
    text(0, L_eixos + 300, 0, 'Y_T', 'FontWeight', 'bold');
    text(0, 0, L_eixos + 300, 'Z_T / Z_I (Polo N)', 'FontWeight', 'bold');

    % Eixos SGI (Inerciais - Cor Vermelha)
    XI_dir = Rz' * [L_eixos; 0; 0];
    YI_dir = Rz' * [0; L_eixos; 0];

    plot3([0 XI_dir(1)], [0 XI_dir(2)], [0 XI_dir(3)], 'r-', 'LineWidth', 1.5);
    plot3([0 YI_dir(1)], [0 YI_dir(2)], [0 YI_dir(3)], 'r-', 'LineWidth', 1.5);
    text(XI_dir(1) + 300, XI_dir(2), XI_dir(3), 'X_I (Aries \Upsilon)', 'Color', 'r', 'FontWeight', 'bold');
    text(YI_dir(1), YI_dir(2) + 300, YI_dir(3), 'Y_I', 'Color', 'r', 'FontWeight', 'bold');

    % Plotagem do Arco do Tempo Sideral de Greenwich (TSG)
    tsg_rad = deg2rad(tsg_graus);
    theta_arc = linspace(0, tsg_rad, 40);
    r_arc = 3500;
    plot3(r_arc * cos(theta_arc), r_arc * sin(theta_arc), zeros(1, 40), 'm-', 'LineWidth', 2);
    text(r_arc * cos(tsg_rad/2) + 200, r_arc * sin(tsg_rad/2) + 200, 300, ...
         sprintf('\\theta_g = %.2f°', tsg_graus), 'Color', 'm', 'FontWeight', 'bold');

    % Ponto da Estação Terrestre
    plot3(aT(1), aT(2), aT(3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(aT(1) + 200, aT(2) + 200, aT(3) + 200, 'Xichang (a_T)', 'Color', 'b', 'FontWeight', 'bold');

    % Vetor Posição da Estação
    plot3([0 aT(1)], [0 aT(2)], [0 aT(3)], 'b-', 'LineWidth', 2);

    % Projeção da Estação para visualização geométrica
    plot3([0 aT(1)], [0 aT(2)], [0 0], 'b--', 'LineWidth', 1);
    plot3([aT(1) aT(1)], [aT(2) aT(2)], [0 aT(3)], 'b:', 'LineWidth', 1.2);

    % Vetor Direção Inercial (Versor aI)
    L_versor = 2000;
    plot3([aT(1), aT(1) + L_versor*aI(1)], ...
          [aT(2), aT(2) + L_versor*aI(2)], ...
          [aT(3), aT(3) + L_versor*aI(3)], ...
          'g-', 'LineWidth', 2.5);
    text(aT(1) + L_versor*aI(1) + 200, ...
         aT(2) + L_versor*aI(2) + 200, ...
         aT(3) + L_versor*aI(3) + 200, ...
         'Zenith / a_I', 'Color', [0, 0.5, 0], 'FontWeight', 'bold');

    title('Geometria do Problema 3: SCT, SGI e Posição da Estação');
    xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
    legend('Terra (Esfera)', 'Equador', 'Meridiano Greenwich', ...
           'SCT Eixo X_T', 'SCT Eixo Y_T', 'SCT Eixo Z_T', ...
           'SGI Eixo X_I', 'SGI Eixo Y_I', 'Arco TSG (\theta_g)', ...
           'Estação (XSLC)', 'Vetor Estação', 'Zenith (a_I)', ...
           'Location', 'northeastoutside');

    hold off;

    % Salvar a imagem gerada
    try
        saveas(gcf, 'geometria_3d_pb3.png');
        saveas(gcf, 'geometria_3d_pb3.svg');
        fprintf('Gráfico salvo como geometria_3d_pb3.png e geometria_3d_pb3.svg no diretório corrente.\n');
    catch
        fprintf('Aviso: Não foi possível salvar o gráfico automaticamente.\n');
    end

    % =========================================================================
    % 8. SEGUNDA VISUALIZAÇÃO GRÁFICA 3D (RENDERIZADA COM TEXTURA DA TERRA)
    % =========================================================================
    if exist('earth_texture.jpg', 'file')
        figure(2);
        hold on; grid on; axis equal; view(135, 25);

        % Esfera com textura fotográfica da Terra (mar e continentes)
        [X_esf, Y_esf, Z_esf] = sphere(50);
        img_earth = imread('earth_texture.jpg');
        img_earth = flipud(img_earth);
        % Alinha a longitude 0° (centro da imagem) com a coluna 1 da esfera (que inicia em 0°)
        img_earth = circshift(img_earth, [0, size(img_earth, 2)/2, 0]);
        surf(X_esf * ae, Y_esf * ae, Z_esf * ae, 'FaceColor', 'texturemap', ...
             'CData', img_earth, 'EdgeColor', 'none', 'FaceAlpha', 0.95);

        % Plotagem do Equador e do Meridiano de Greenwich
        h_eq = plot3(ae*cos(t_circle), ae*sin(t_circle), zeros(1, 100), 'k--', 'LineWidth', 1.0); % Equador
        h_gr = plot3(ae*cos(t_circle), zeros(1, 100), ae*sin(t_circle), 'Color', [0.3, 0.3, 0.3], 'LineStyle', ':', 'LineWidth', 1.0); % Greenwich

        % Eixos SCT (Fixos na Terra - Cor Preta)
        h_sct = plot3([0 L_eixos], [0 0], [0 0], 'k-', 'LineWidth', 1.5);
        plot3([0 0], [0 L_eixos], [0 0], 'k-', 'LineWidth', 1.5);
        plot3([0 0], [0 0], [0 L_eixos], 'k-', 'LineWidth', 1.5);
        text(L_eixos + 300, 0, 0, 'X_T (Greenwich)', 'FontWeight', 'bold');
        text(0, L_eixos + 300, 0, 'Y_T', 'FontWeight', 'bold');
        text(0, 0, L_eixos + 300, 'Z_T / Z_I (Polo N)', 'FontWeight', 'bold');

        % Eixos SGI (Inerciais - Cor Vermelha)
        h_sgi = plot3([0 XI_dir(1)], [0 XI_dir(2)], [0 XI_dir(3)], 'r-', 'LineWidth', 1.5);
        plot3([0 YI_dir(1)], [0 YI_dir(2)], [0 YI_dir(3)], 'r-', 'LineWidth', 1.5);
        text(XI_dir(1) + 300, XI_dir(2), XI_dir(3), 'X_I (Aries \Upsilon)', 'Color', 'r', 'FontWeight', 'bold');
        text(YI_dir(1), YI_dir(2) + 300, YI_dir(3), 'Y_I', 'Color', 'r', 'FontWeight', 'bold');

        % Arco do TSG
        h_arc = plot3(r_arc * cos(theta_arc), r_arc * sin(theta_arc), zeros(1, 40), 'm-', 'LineWidth', 2);
        text(r_arc * cos(tsg_rad/2) + 200, r_arc * sin(tsg_rad/2) + 200, 300, ...
             sprintf('\\theta_g = %.2f°', tsg_graus), 'Color', 'm', 'FontWeight', 'bold');

        % Ponto da Estação Terrestre (XSLC)
        h_est = plot3(aT(1), aT(2), aT(3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        text(aT(1) + 200, aT(2) + 200, aT(3) + 200, 'Xichang (a_T)', 'Color', 'b', 'FontWeight', 'bold');

        % Vetores de Posição e Projeções
        h_vest = plot3([0 aT(1)], [0 aT(2)], [0 aT(3)], 'b-', 'LineWidth', 2);
        plot3([0 aT(1)], [0 aT(2)], [0 0], 'b--', 'LineWidth', 1);
        plot3([aT(1) aT(1)], [aT(2) aT(2)], [0 aT(3)], 'b:', 'LineWidth', 1.2);

        % Vetor Zenital Inercial (aI)
        h_zen = plot3([aT(1), aT(1) + L_versor*aI(1)], ...
              [aT(2), aT(2) + L_versor*aI(2)], ...
              [aT(3), aT(3) + L_versor*aI(3)], ...
              'g-', 'LineWidth', 2.5);
        text(aT(1) + L_versor*aI(1) + 200, ...
             aT(2) + L_versor*aI(2) + 200, ...
             aT(3) + L_versor*aI(3) + 200, ...
             'Zenith / a_I', 'Color', [0, 0.5, 0], 'FontWeight', 'bold');

        title('Geometria do Problema 3: SCT, SGI e Posição da Estação (Renderizado)');
        xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
        legend([h_eq, h_gr, h_sct, h_sgi, h_arc, h_est, h_vest, h_zen], ...
               'Equador', 'Meridiano Greenwich', 'SCT Eixos', 'SGI Eixos', ...
               'Arco TSG (\theta_g)', 'Estação (XSLC)', 'Vetor Estação', 'Zenith (a_I)', ...
               'Location', 'northeastoutside');

        hold off;

        % Salvar a imagem renderizada
        try
            saveas(gcf, 'geometria_3d_pb3_renderizado.png');
            saveas(gcf, 'geometria_3d_pb3_renderizado.svg');
            fprintf('Gráfico renderizado salvo como geometria_3d_pb3_renderizado.png e geometria_3d_pb3_renderizado.svg no diretório corrente.\n');
        catch
            fprintf('Aviso: Não foi possível salvar o gráfico renderizado automaticamente.\n');
        end
    end
end

% =========================================================================
% FUNÇÕES AUXILIARES LOCAIS
% =========================================================================

function [JD, DJM0, DJM_frac] = DiaJulianoLocal(ano, mes, dia, hora, minuto, segundo)
    if nargin < 4, hora = 0; end
    if nargin < 5, minuto = 0; end
    if nargin < 6, segundo = 0; end

    I = ano;
    J = mes;
    K = dia;

    dj1 = K - 32075 + fix((1461 * (I + 4800 + fix((J - 14)/12))) / 4);
    dj2 = fix((367 * (J - 2 - 12 * fix((J - 14)/12))) / 12);
    dj3 = fix((-3 * fix((I + 4900 + fix((J - 14)/12)) / 100)) / 4);

    JD = dj1 + dj2 + dj3;
    DJM0 = JD - 2400001;
    frac_dia = (hora + minuto/60 + segundo/3600) / 24;
    DJM_frac = DJM0 + frac_dia;
end

function tsg = TSGLocal(ano, mes, dia, hora, minuto, segundo)
    DJ0 = DiaJulianoLocal(ano, mes, dia, 0, 0, 0) - 0.5;
    delta_t = hora * 3600 + minuto * 60 + segundo;

    Tu = (DJ0 - 2415020.0) / 36525;
    L = 365.24219879 - 6.14e-6 * Tu;
    d_theta_dt = (1/240) * (1 + 1/L);
    theta_g0 = 99.69098329 + 36000.76893 * Tu + 3.87080e-4 * Tu^2;
    tsg = mod(theta_g0 + d_theta_dt * delta_t, 360);
end
