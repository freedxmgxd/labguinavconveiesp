% Script para testar DiaJuliano e TSG
clear; clc;

disp('=== TESTE: DATA JULIANA ===')

% 1) Dias entre 12/Dez/2008 e 25/Mar/2020
jd_ini = DiaJuliano(2008, 12, 12);
jd_fim = DiaJuliano(2020, 3, 25);
fprintf('Pergunta 1: %d - %d = %d dias\n\n', jd_fim, jd_ini, jd_fim - jd_ini)

% 2) Exemplo 1 do PDF: 25/Dez/1981
[jd1, djm0_1, ~] = DiaJuliano(1981, 12, 25);
fprintf('Exemplo 1 (25/12/1981):\n')
fprintf('  JD  = %d  (esperado: 2444964)\n', jd1)
fprintf('  DJM0= %d  (esperado: 44963)\n\n', djm0_1)

% 3) Com fração de hora: 09:10:35
[~, ~, djm_frac1] = DiaJuliano(1981, 12, 25, 9, 10, 35);
fprintf('Exemplo 1 com 09h10m35s:\n')
fprintf('  DJM_frac = %.15f\n', djm_frac1)
fprintf('  esperado = 44963.382349537037\n\n')

% 4) Exemplo 2: 24/Mar/2020
[jd2, djm0_2, ~] = DiaJuliano(2020, 3, 24);
fprintf('Exemplo 2 (24/03/2020):\n')
fprintf('  JD  = %d  (esperado: 2458933)\n', jd2)
fprintf('  DJM0= %d  (esperado: 58932)\n\n', djm0_2)

% 5) Data da aula: 01/Jul/2026 14:11:11
[jd_aula, djm0_aula, djm_frac_aula] = DiaJuliano(2026, 7, 1, 14, 11, 11);
fprintf('Aula (01/07/2026 14:11:11):\n')
fprintf('  JD       = %d\n', jd_aula)
fprintf('  DJM0     = %d\n', djm0_aula)
fprintf('  DJM_frac = %.15f\n\n', djm_frac_aula)

disp('=== TESTE: TSG ===')

% Exemplo do PDF (Seção 5): 1962/Out/12, 10h15m30s UT
tsg_ex = TSG(1962, 10, 12, 10, 15, 30);
fprintf('Exemplo PDF (12/10/1962 10h15m30s):\n')
fprintf('  TSG = %.4f graus  (esperado: 174.3880)\n', tsg_ex)

% Tempo Sideral Local em Wantig (lambda_E = 298.2213 graus)
tsl = mod(tsg_ex + 298.2213, 360);
fprintf('  TSL = %.4f graus  (esperado: 112.6093)\n\n', tsl)

% TSG da aula
tsg_aula = TSG(2026, 7, 1, 14, 11, 11);
fprintf('Aula (01/07/2026 14:11:11):\n')
fprintf('  TSG = %.4f graus\n', tsg_aula)
