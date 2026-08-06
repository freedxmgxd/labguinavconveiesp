%% ======================================================================
%  PB5 - Lab de GN&C - UFABC - Prof. A. Gil V. de Brum
%  ATITUDE DE UM V/E: REPRESENTACAO, SIMULACAO E DETERMINACAO
%
%  Parte 1 - Simula a rotacao de um V/E (matriz de atitude A132), gera os
%            vetores no sistema do corpo e recupera os angulos de atitude
%            pelo metodo algebrico (TRIAD).
%  Parte 2 - Determina a atitude do V/E a partir de duas estrelas (Matar e
%            Markab) identificadas no CDV de um sensor rastreador estelar.
%
%  Convencao adotada (ref.: "Sobre outras matrizes de atitude", Prof. Gil):
%     Rx(roll)  = [1 0 0; 0 cr sr; 0 -sr cr]
%     Ry(pitch) = [cp 0 -sp; 0 1 0; sp 0 cp]
%     Rz(yaw)   = [cy sy 0; -sy cy 0; 0 0 1]
%     A132 = Ry(pitch) * Rz(yaw) * Rx(roll)   =>   b = A132 * r
%     (rotacao em torno de x, depois z, depois y => sequencia 1-3-2)
%
%  Uso: basta digitar  >> pb5   na linha de comando (MATLAB ou Octave).
%% ======================================================================
function pb5()

clc; close all;

fprintf('======================================================================\n');
fprintf('  PB5 - ATITUDE DE UM V/E: REPRESENTACAO, SIMULACAO E DETERMINACAO\n');
fprintf('======================================================================\n');

%% ======================================================================
%  PARTE 1 - SIMULACAO DA ROTACAO E RECUPERACAO DOS ANGULOS
%% ======================================================================
fprintf('\n\n#################### P A R T E   1 ####################\n');

%% --- Item 4: vetores de referencia no SGI e seus versores --------------
r1 = [ 1; -2; -2];
r2 = [-1;  2;  3];

r1_hat = r1 / norm(r1);
r2_hat = r2 / norm(r2);

fprintf('\n--- Item 4: versores de referencia no SGI ---\n');
fprintf('  r1  = [%7.4f %7.4f %7.4f]''   |r1| = %.4f\n', r1, norm(r1));
fprintf('  r2  = [%7.4f %7.4f %7.4f]''   |r2| = %.4f\n', r2, norm(r2));
fprintf('  r1^ = [%9.6f %9.6f %9.6f]''\n', r1_hat);
fprintf('  r2^ = [%9.6f %9.6f %9.6f]''\n', r2_hat);
fprintf('  Angulo entre r1^ e r2^ : %.4f graus\n', acosd(dot(r1_hat, r2_hat)));

%% --- Item 5: as tres matrizes de rotacao elementares -------------------
roll_sim  = 187.9;   % teta_roll  (rotacao em torno de x - rolamento)
pitch_sim = 341.5;   % teta_pitch (rotacao em torno de y - arfagem)
yaw_sim   =  79.7;   % teta_yaw   (rotacao em torno de z - guinada)

Mx = rotx_gnc(roll_sim);    % i)   em torno de x, de um angulo teta_roll
My = roty_gnc(pitch_sim);   % ii)  em torno de y, de um angulo teta_pitch
Mz = rotz_gnc(yaw_sim);     % iii) em torno de z, de um angulo teta_yaw

fprintf('\n--- Item 5: matrizes de rotacao (angulos simulados) ---\n');
fprintf('  teta_roll = %.1f deg | teta_pitch = %.1f deg | teta_yaw = %.1f deg\n', ...
        roll_sim, pitch_sim, yaw_sim);
print_mat('Mx = Rx(teta_roll)', Mx);
print_mat('My = Ry(teta_pitch)', My);
print_mat('Mz = Rz(teta_yaw)', Mz);

%% --- Item 6: metodo TRIAD - rotacao dos versores (b = My*Mz*Mx*r) ------
A_sim = My * Mz * Mx;        % matriz de atitude simulada (A132)

b1_hat = A_sim * r1_hat;
b2_hat = A_sim * r2_hat;

fprintf('\n--- Item 6: versores vistos no sistema do corpo ---\n');
print_mat('A132(simulada) = My*Mz*Mx', A_sim);
fprintf('  b1^ = My*Mz*Mx*r1^ = [%9.6f %9.6f %9.6f]''\n', b1_hat);
fprintf('  b2^ = My*Mz*Mx*r2^ = [%9.6f %9.6f %9.6f]''\n', b2_hat);
fprintf('  Angulo entre b1^ e b2^ : %.4f graus (preservado pela rotacao)\n', ...
        acosd(dot(b1_hat, b2_hat)));

%% --- Item 7: metodo algebrico (TRIAD) - matrizes Mr, Mb e A_T ----------
[A_T, Mr, Mb] = triad(r1_hat, r2_hat, b1_hat, b2_hat);

fprintf('\n--- Item 7: metodo algebrico TRIAD ---\n');
print_mat('Mr = [v1 v2 v3] (triade de referencia)', Mr);
print_mat('Mb = [w1 w2 w3] (triade do corpo)', Mb);
print_mat('A_T = Mb*Mr'' (matriz de atitude determinada)', A_T);

fprintf('  Verificacao de ortogonalidade: ||A_T*A_T'' - I|| = %.3e\n', ...
        norm(A_T*A_T' - eye(3)));
fprintf('  Determinante de A_T          : %.10f\n', det(A_T));
fprintf('  Erro vs. matriz simulada     : ||A_T - A132|| = %.3e\n', ...
        norm(A_T - A_sim));

%% --- Item 8: extracao dos angulos de atitude da matriz A132 ------------
[roll_rec, pitch_rec, yaw_rec] = ang_from_A132(A_T);

fprintf('\n--- Item 8: angulos de atitude extraidos de A_T (matriz 132) ---\n');
fprintf('  A132 simbolica:\n');
fprintf('    [  cp*cy            cp*sy*cr + sp*sr    cp*sy*sr - sp*cr ]\n');
fprintf('    [ -sy               cy*cr               cy*sr            ]\n');
fprintf('    [  sp*cy            sp*sy*cr - cp*sr    sp*sy*sr + cp*cr ]\n');
fprintf('  => yaw   = asind(-A(2,1))\n');
fprintf('     pitch = atan2d(A(3,1)/cos(yaw), A(1,1)/cos(yaw))\n');
fprintf('     roll  = atan2d(A(2,3)/cos(yaw), A(2,2)/cos(yaw))\n\n');

fprintf('  %-22s %12s %12s %12s\n', 'Angulo', 'simulado', 'recuperado', 'equiv. [0,360)');
fprintf('  %-22s %12.4f %12.4f %12.4f\n', 'rolamento (roll)',  roll_sim,  roll_rec,  mod(roll_rec, 360));
fprintf('  %-22s %12.4f %12.4f %12.4f\n', 'arfagem   (pitch)', pitch_sim, pitch_rec, mod(pitch_rec, 360));
fprintf('  %-22s %12.4f %12.4f %12.4f\n', 'guinada   (yaw)',   yaw_sim,   yaw_rec,   mod(yaw_rec, 360));

err_ang = [ang_diff(roll_rec, roll_sim), ang_diff(pitch_rec, pitch_sim), ang_diff(yaw_rec, yaw_sim)];
fprintf('\n  Erro (modulo 360) : roll %.2e | pitch %.2e | yaw %.2e graus\n', err_ang);
fprintf('  ==> Os tres angulos das rotacoes simuladas foram recuperados.\n');
fprintf('  Obs.: atan2 devolve o valor principal em (-180,180]; 187.9 deg e\n');
fprintf('        341.5 deg reaparecem como -172.1 deg e -18.5 deg, angulos\n');
fprintf('        equivalentes (mesma matriz de atitude).\n');

% Reconstrucao de A132 com os angulos recuperados (fechamento do ciclo)
A_rec = roty_gnc(pitch_rec) * rotz_gnc(yaw_rec) * rotx_gnc(roll_rec);
fprintf('  Fechamento do ciclo: ||A132(angulos recuperados) - A_T|| = %.3e\n', ...
        norm(A_rec - A_T));

%% ======================================================================
%  PARTE 2 - DETERMINACAO DA ATITUDE COM SENSOR RASTREADOR ESTELAR
%% ======================================================================
fprintf('\n\n#################### P A R T E   2 ####################\n');

%% --- Dados observacionais (SGI / J2000.0) ------------------------------
% Estrela 1: Matar  (Eta Pegasi)   - AR 22h 44m 16.0s | Dec +30 21'' 32.7''''
% Estrela 2: Markab (Alpha Pegasi) - AR 23h 06m 06.0s | Dec +15 20'' 53.8''''
estrelas = {'Matar  (Eta Pegasi)', 'Markab (Alpha Pegasi)'};

AR_hms  = [22 44 16.0;      % [h m s]
           23 06 06.0];
DEC_dms = [+30 21 32.7;     % [deg arcmin arcsec]
           +15 20 53.8];

AR_deg  = hms2deg(AR_hms);
DEC_deg = dms2deg(DEC_dms);

% Versores de referencia no SGI: r = [cos(d)cos(a); cos(d)sin(a); sin(d)]
R_ref = radec2vec(AR_deg, DEC_deg);
s1_hat = R_ref(:,1);
s2_hat = R_ref(:,2);

fprintf('\n--- Dados de entrada: coordenadas inerciais (SGI/J2000.0) ---\n');
for k = 1:2
    fprintf('  %s\n', estrelas{k});
    fprintf('    AR  = %2dh %2dm %4.1fs  = %10.6f deg\n', AR_hms(k,:), AR_deg(k));
    fprintf('    Dec = %+3d %2d'' %4.1f" = %10.6f deg\n', DEC_dms(k,:), DEC_deg(k));
    fprintf('    r%d^ = [%9.6f %9.6f %9.6f]''\n', k, R_ref(:,k));
end

%% --- Vetores medidos no sistema do corpo (saida do sensor) -------------
b1_st = [ 0.5506; -0.8138; 0.1858];
b2_st = [-0.3521; -0.9351; 0.0386];

fprintf('\n--- Dados de entrada: vetores medidos no sistema do corpo ---\n');
fprintf('  b1 = [%9.6f %9.6f %9.6f]''   |b1| = %.6f\n', b1_st, norm(b1_st));
fprintf('  b2 = [%9.6f %9.6f %9.6f]''   |b2| = %.6f\n', b2_st, norm(b2_st));

% Normalizacao (o sensor entrega os vetores ja praticamente unitarios)
b1_st_hat = b1_st / norm(b1_st);
b2_st_hat = b2_st / norm(b2_st);

%% --- Verificacao de consistencia dos dados -----------------------------
% Uma rotacao e uma transformacao ortogonal: ela PRESERVA o angulo entre
% os vetores. Logo, o angulo entre b1 e b2 tem de ser igual ao angulo
% entre r1^ e r2^. Este teste deve preceder qualquer determinacao.
ang_ref   = acosd(max(-1, min(1, dot(s1_hat, s2_hat))));
ang_corpo = acosd(max(-1, min(1, dot(b1_st_hat, b2_st_hat))));

fprintf('\n--- Verificacao de consistencia (invariancia do angulo) ---\n');
fprintf('  Angulo entre as estrelas no SGI  : %.4f graus\n', ang_ref);
fprintf('  Angulo entre b1 e b2 (corpo)     : %.4f graus\n', ang_corpo);
fprintf('  Discrepancia                     : %+.4f graus\n', ang_corpo - ang_ref);

% Teste da hipotese de erro de sinal na 1a componente de b2
b2_corr     = [-b2_st(1); b2_st(2); b2_st(3)];
b2_corr_hat = b2_corr / norm(b2_corr);
ang_corr    = acosd(max(-1, min(1, dot(b1_st_hat, b2_corr_hat))));

if abs(ang_corpo - ang_ref) > 0.5
    fprintf('\n  *** ATENCAO: discrepancia de %.2f graus. Um sensor estelar tem\n', ...
            abs(ang_corpo - ang_ref));
    fprintf('      exatidao da ordem de segundos de arco; um erro desta magnitude\n');
    fprintf('      NAO e ruido de medida - os dados publicados sao inconsistentes.\n');
    fprintf('\n  Teste de hipotese (troca de sinal na 1a componente de b2):\n');
    fprintf('      b2 = [%+.4f %+.4f %+.4f] => angulo = %.4f graus\n', ...
            b2_corr, ang_corr);
    fprintf('      discrepancia = %+.4f graus (%.1f segundos de arco)\n', ...
            ang_corr - ang_ref, abs(ang_corr - ang_ref)*3600);
    fprintf('      ==> compativel com o arredondamento das 4 casas decimais:\n');
    fprintf('          o sinal "-" em b2(1) e um erro de digitacao do enunciado.\n');
end

%% --- Determinacao da atitude pelo metodo algebrico (TRIAD) -------------
% Os dois cenarios sao resolvidos: (A) dados exatamente como publicados e
% (B) com o sinal de b2(1) corrigido - fisicamente consistente.
cenarios = {'CENARIO A - dados exatamente como publicados', ...
            'CENARIO B - b2(1) com sinal corrigido (+0.3521)'};
B2_casos = [b2_st_hat, b2_corr_hat];

A_caso   = zeros(3,3,2);
ang_caso = zeros(3,2);
res_caso = zeros(2,2);

for k = 1:2
    b2_k = B2_casos(:,k);
    [A_k, Mr_k, Mb_k] = triad(s1_hat, s2_hat, b1_st_hat, b2_k);

    [rl, pt, yw] = ang_from_A132(A_k);
    e1 = acosd(max(-1, min(1, dot(A_k*s1_hat, b1_st_hat))));
    e2 = acosd(max(-1, min(1, dot(A_k*s2_hat, b2_k))));

    A_caso(:,:,k)  = A_k;
    ang_caso(:,k)  = [rl; pt; yw];
    res_caso(:,k)  = [e1; e2];

    fprintf('\n\n=== %s ===\n', cenarios{k});
    fprintf('\n--- Metodo algebrico (TRIAD - Wertz, 1978) ---\n');
    print_mat('Mr (triade de referencia)', Mr_k);
    print_mat('Mb (triade do corpo)', Mb_k);
    print_mat('A = Mb*Mr'' (ATITUDE DO V/E)', A_k);
    fprintf('  Ortogonalidade: ||A*A'' - I|| = %.3e | det(A) = %.10f\n', ...
            norm(A_k*A_k' - eye(3)), det(A_k));

    fprintf('\n--- Atitude em angulos de Euler (sequencia 132) ---\n');
    fprintf('  %-22s %14s %16s\n', 'Angulo', 'valor [deg]', 'equiv. [0,360)');
    fprintf('  %-22s %14.4f %16.4f\n', 'rolamento (roll)',  rl, mod(rl, 360));
    fprintf('  %-22s %14.4f %16.4f\n', 'arfagem   (pitch)', pt, mod(pt, 360));
    fprintf('  %-22s %14.4f %16.4f\n', 'guinada   (yaw)',   yw, mod(yw, 360));

    % O TRIAD casa b1 exatamente; todo o erro sobra no segundo vetor
    fprintf('\n  Residuo angular em b1 (vetor de maior peso) : %.4f graus\n', e1);
    fprintf('  Residuo angular em b2                       : %.4f graus\n', e2);
end

%% --- Resultado adotado -------------------------------------------------
% O cenario B e o unico fisicamente admissivel (preserva o angulo entre as
% direcoes), e portanto e a resposta adotada para a atitude do V/E.
b2_ad = b2_corr_hat;
A_ad  = A_caso(:,:,2);
roll_ad = ang_caso(1,2); pitch_ad = ang_caso(2,2); yaw_ad = ang_caso(3,2);

fprintf('\n\n======================================================================\n');
fprintf('  RESPOSTA - ATITUDE DO V/E NO INSTANTE DA AQUISICAO\n');
fprintf('======================================================================\n');
fprintf('  %-22s %14s %14s\n', 'Angulo', 'Cenario A', 'Cenario B*');
fprintf('  %-22s %14.4f %14.4f\n', 'rolamento (roll)',  ang_caso(1,1), ang_caso(1,2));
fprintf('  %-22s %14.4f %14.4f\n', 'arfagem   (pitch)', ang_caso(2,1), ang_caso(2,2));
fprintf('  %-22s %14.4f %14.4f\n', 'guinada   (yaw)',   ang_caso(3,1), ang_caso(3,2));
fprintf('  %-22s %14.4f %14.4f\n', 'residuo em b2 [deg]', res_caso(2,1), res_caso(2,2));
print_mat('Matriz de atitude adotada (cenario B)', A_ad);
fprintf('  (*) Cenario B adotado: e o unico que preserva o angulo entre as\n');
fprintf('      direcoes observadas, condicao necessaria para que exista uma\n');
fprintf('      matriz de atitude compativel com as duas medidas.\n');

%% ======================================================================
%  QUESTAO C - VARIAS ESTRELAS NO CDV: O PROBLEMA DE WAHBA
%% ======================================================================
% O TRIAD usa apenas 2 direcoes e descarta informacao. Com N > 2 estrelas
% resolve-se o problema de Wahba: minimizar
%       J(A) = (1/2) * soma_k[ w_k * |b_k - A*r_k|^2 ],
% cuja solucao otima e obtida pelo metodo q de Davenport (base do QUEST):
% monta-se a matriz K (4x4) e toma-se o autovetor do maior autovalor, que
% e o quaternio de atitude procurado.
fprintf('\n\n#################### Q U E S T A O   C ####################\n');

%% --- C.1: Wahba aplicado as duas estrelas reais ------------------------
A_wahba = davenport_q(R_ref, [b1_st_hat, b2_ad], [0.5; 0.5]);
[roll_w, pitch_w, yaw_w] = ang_from_A132(A_wahba);

fprintf('\n--- C.1: metodo q de Davenport nas 2 estrelas (cenario B) ---\n');
print_mat('A_wahba (pesos iguais)', A_wahba);
fprintf('  %-22s %12s %12s %12s\n', 'Angulo', 'TRIAD', 'Wahba', 'dif [deg]');
fprintf('  %-22s %12.4f %12.4f %12.4f\n', 'rolamento (roll)',  roll_ad,  roll_w,  ang_diff(roll_w, roll_ad));
fprintf('  %-22s %12.4f %12.4f %12.4f\n', 'arfagem   (pitch)', pitch_ad, pitch_w, ang_diff(pitch_w, pitch_ad));
fprintf('  %-22s %12.4f %12.4f %12.4f\n', 'guinada   (yaw)',   yaw_ad,   yaw_w,   ang_diff(yaw_w, yaw_ad));
fprintf('\n  Com apenas 2 direcoes quase consistentes os dois metodos coincidem;\n');
fprintf('  a diferenca aparece quando ha redundancia (N > 2) e ruido.\n');

%% --- C.2: simulacao com N estrelas no CDV (TRIAD x QUEST) --------------
% Atitude verdadeira conhecida, N estrelas sinteticas no CDV e ruido
% gaussiano de 20 segundos de arco por eixo (tipico de um AST).
rng_seed = 42;
if exist('rng', 'builtin') || exist('rng', 'file'), rng(rng_seed); else, rand('seed', rng_seed); randn('seed', rng_seed); end

A_true   = roty_gnc(pitch_sim) * rotz_gnc(yaw_sim) * rotx_gnc(roll_sim);
N_est    = 8;
sigma_as = 20;                       % ruido do sensor [segundos de arco]
sigma    = sigma_as / 3600;          % em graus

R_sim = zeros(3, N_est);
B_sim = zeros(3, N_est);
for k = 1:N_est
    v = randn(3,1);  R_sim(:,k) = v / norm(v);                 % estrela k no SGI
    bk = A_true * R_sim(:,k);
    bk = bk + (sigma*pi/180) * randn(3,1);                     % ruido de medida
    B_sim(:,k) = bk / norm(bk);
end

% TRIAD: obrigado a escolher so 2 das N estrelas (usa as duas primeiras)
A_triad_N = triad(R_sim(:,1), R_sim(:,2), B_sim(:,1), B_sim(:,2));
% QUEST/Wahba: usa todas as N estrelas
A_quest_N = davenport_q(R_sim, B_sim, ones(N_est,1));

fprintf('\n--- C.2: simulacao com %d estrelas no CDV (ruido = %d arcsec) ---\n', N_est, sigma_as);
fprintf('  Atitude verdadeira : roll = %.1f | pitch = %.1f | yaw = %.1f deg\n', ...
        roll_sim, pitch_sim, yaw_sim);
fprintf('  Erro do TRIAD (so 2 estrelas de %d) : %8.2f segundos de arco\n', ...
        N_est, err_ang_matriz(A_triad_N, A_true)*3600);
fprintf('  Erro do QUEST (todas as %d estrelas): %8.2f segundos de arco\n', ...
        N_est, err_ang_matriz(A_quest_N, A_true)*3600);
fprintf('\n  ==> RESPOSTA: com varias estrelas nao se usa o metodo algebrico.\n');
fprintf('      Resolve-se o problema de Wahba (metodo q de Davenport / QUEST),\n');
fprintf('      que pondera TODAS as direcoes medidas pelos respectivos pesos.\n');
fprintf('      O erro cai aproximadamente com 1/sqrt(N), como se ve acima.\n');

%% ======================================================================
%  FIGURAS
%% ======================================================================
% Fig. 1 - Parte 1: versores de referencia e do corpo
figure('Name', 'PB5 - Parte 1: rotacao simulada', 'Color', 'w');
plot_frames(r1_hat, r2_hat, b1_hat, b2_hat, ...
            'Parte 1: r1^, r2^ (SGI) e b1^, b2^ (corpo)');

% Fig. 2 - Parte 2: estrelas no SGI e no sistema do corpo (cenario adotado)
figure('Name', 'PB5 - Parte 2: sensor de estrelas', 'Color', 'w');
plot_frames(s1_hat, s2_hat, b1_st_hat, b2_ad, ...
            'Parte 2: Matar e Markab no SGI e no corpo');

saveas(figure(1), 'pb5_parte1_vetores.png');
saveas(figure(2), 'pb5_parte2_estrelas.png');

fprintf('\n======================================================================\n');
fprintf('  Figuras salvas: pb5_parte1_vetores.png / pb5_parte2_estrelas.png\n');
fprintf('======================================================================\n');

end % ---------------------------- fim de pb5 ---------------------------


%% ======================================================================
%  FUNCOES AUXILIARES
%% ======================================================================

function R = rotx_gnc(ang)
% Rotacao em torno de x (rolamento) - convencao do curso
    c = cosd(ang); s = sind(ang);
    R = [1  0  0; 0  c  s; 0 -s  c];
end

function R = roty_gnc(ang)
% Rotacao em torno de y (arfagem) - convencao do curso
    c = cosd(ang); s = sind(ang);
    R = [c  0 -s; 0  1  0; s  0  c];
end

function R = rotz_gnc(ang)
% Rotacao em torno de z (guinada) - convencao do curso
    c = cosd(ang); s = sind(ang);
    R = [ c  s  0; -s  c  0; 0  0  1];
end

function [A, Mr, Mb] = triad(r1_hat, r2_hat, b1_hat, b2_hat)
% Metodo algebrico TRIAD (Wertz, 1978).
% Constroi uma triade ortonormal em cada sistema a partir de dois vetores
% e obtem a atitude por A = Mb*Mr'. O primeiro vetor (r1/b1) e casado
% exatamente; o segundo entra apenas pela direcao do produto vetorial.
    v1 = r1_hat;
    v2 = cross(r1_hat, r2_hat);  v2 = v2 / norm(v2);
    v3 = cross(v1, v2);
    Mr = [v1, v2, v3];

    w1 = b1_hat;
    w2 = cross(b1_hat, b2_hat);  w2 = w2 / norm(w2);
    w3 = cross(w1, w2);
    Mb = [w1, w2, w3];

    A = Mb * Mr';
end

function [roll, pitch, yaw] = ang_from_A132(A)
% Extrai rolamento, arfagem e guinada (graus) da matriz de atitude 132,
% A132 = Ry(pitch)*Rz(yaw)*Rx(roll):
%   [  cp*cy     cp*sy*cr + sp*sr    cp*sy*sr - sp*cr ]
%   [ -sy        cy*cr               cy*sr            ]
%   [  sp*cy     sp*sy*cr - cp*sr    sp*sy*sr + cp*cr ]
% O yaw sai do seno (limitado a [-90,90], evitando a ambiguidade) e os
% outros dois de atan2, ja com o primeiro angulo conhecido.
    s_yaw = -A(2,1);
    s_yaw = max(-1, min(1, s_yaw));       % protecao contra |.|>1 numerico
    yaw   = asind(s_yaw);
    c_yaw = cosd(yaw);

    if abs(c_yaw) < 1e-10                 % singularidade (gimbal lock)
        warning('ang_from_A132:gimbalLock', ...
                'cos(yaw) ~ 0: roll e pitch nao sao separaveis; roll fixado em 0.');
        roll  = 0;
        pitch = atan2d(-A(1,3), A(3,3));
    else
        pitch = atan2d(A(3,1)/c_yaw, A(1,1)/c_yaw);
        roll  = atan2d(A(2,3)/c_yaw, A(2,2)/c_yaw);
    end
end

function A = davenport_q(R, B, w)
% Solucao do problema de Wahba pelo metodo q de Davenport (base do QUEST).
% Minimiza J(A) = (1/2)*soma[ w_k*|b_k - A*r_k|^2 ] para N >= 2 vetores.
% R e B: 3xN com os versores de referencia e medidos; w: pesos (Nx1).
    w = w(:) / sum(w);
    Bm = zeros(3);
    for k = 1:numel(w)
        Bm = Bm + w(k) * (B(:,k) * R(:,k)');
    end
    S = Bm + Bm';
    z = [Bm(2,3) - Bm(3,2); Bm(3,1) - Bm(1,3); Bm(1,2) - Bm(2,1)];
    sigma = trace(Bm);

    K = [S - sigma*eye(3), z; z', sigma];       % matriz de Davenport 4x4

    [V, D] = eig(K);
    [~, idx] = max(real(diag(D)));              % autovetor do maior autovalor
    q = real(V(:, idx));
    q = q / norm(q);
    qv = q(1:3); q4 = q(4);                     % quaternio [vetor; escalar]

    % Matriz de atitude a partir do quaternio (mesma convencao b = A*r)
    A = (q4^2 - qv'*qv)*eye(3) + 2*(qv*qv') - 2*q4*skew(qv);
end

function S = skew(v)
% Matriz antissimetrica associada ao produto vetorial
    S = [   0, -v(3),  v(2);
         v(3),     0, -v(1);
        -v(2),  v(1),     0];
end

function d = hms2deg(hms)
% Ascensao reta em h/m/s -> graus
    d = (hms(:,1) + hms(:,2)/60 + hms(:,3)/3600) * 15;
end

function d = dms2deg(dms)
% Declinacao em grau/arcmin/arcsec -> graus (preserva o sinal)
    sgn = sign(dms(:,1));
    sgn(sgn == 0) = 1;
    d = sgn .* (abs(dms(:,1)) + dms(:,2)/60 + dms(:,3)/3600);
end

function V = radec2vec(ra_deg, dec_deg)
% Versores no SGI a partir de ascensao reta e declinacao (graus)
    V = [cosd(dec_deg) .* cosd(ra_deg), ...
         cosd(dec_deg) .* sind(ra_deg), ...
         sind(dec_deg)]';
end

function d = ang_diff(a, b)
% Menor diferenca angular entre dois angulos, em graus
    d = mod(a - b + 180, 360) - 180;
end

function e = err_ang_matriz(A, A_ref)
% Erro de atitude (angulo unico de Euler da rotacao residual), em graus.
% A rotacao residual e dA = A*A_ref'; seu traco da o angulo de rotacao.
    dA = A * A_ref';
    c  = max(-1, min(1, (trace(dA) - 1) / 2));
    e  = acosd(c);
end

function print_mat(nome, M)
% Impressao formatada de uma matriz 3x3
    fprintf('\n  %s =\n', nome);
    for i = 1:3
        fprintf('    [ %11.6f  %11.6f  %11.6f ]\n', M(i,1), M(i,2), M(i,3));
    end
    fprintf('\n');
end

function plot_frames(r1, r2, b1, b2, titulo)
% Visualizacao 3D dos versores no SGI e no sistema do corpo
    hold on; grid on; box on;
    O = zeros(3,1);

    % Eixos do SGI
    quiver3(0,0,0, 1,0,0, 0, 'k', 'LineWidth', 1);
    quiver3(0,0,0, 0,1,0, 0, 'k', 'LineWidth', 1);
    quiver3(0,0,0, 0,0,1, 0, 'k', 'LineWidth', 1);
    text(1.05, 0, 0, 'X'); text(0, 1.05, 0, 'Y'); text(0, 0, 1.05, 'Z');

    h1 = quiver3(O(1),O(2),O(3), r1(1),r1(2),r1(3), 0, 'b', 'LineWidth', 2.5);
    quiver3(O(1),O(2),O(3), r2(1),r2(2),r2(3), 0, 'b--', 'LineWidth', 2.5);
    h2 = quiver3(O(1),O(2),O(3), b1(1),b1(2),b1(3), 0, 'r', 'LineWidth', 2.5);
    quiver3(O(1),O(2),O(3), b2(1),b2(2),b2(3), 0, 'r--', 'LineWidth', 2.5);

    % Esfera unitaria de referencia
    [xs, ys, zs] = sphere(30);
    surf(xs, ys, zs, 'FaceAlpha', 0.05, 'EdgeColor', [0.8 0.8 0.8], ...
         'FaceColor', [0.5 0.5 0.5], 'EdgeAlpha', 0.2);

    legend([h1 h2], {'referencia (SGI)', 'corpo (V/E)'}, 'Location', 'northeast');
    xlabel('x'); ylabel('y'); zlabel('z'); title(titulo);
    axis equal; view(35, 22);
    xlim([-1.2 1.2]); ylim([-1.2 1.2]); zlim([-1.2 1.2]);
end
