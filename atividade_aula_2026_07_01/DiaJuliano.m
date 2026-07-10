function [JD, DJM0, DJM_frac] = DiaJuliano(ano, mes, dia, hora, minuto, segundo)
% Função para calcular a Data Juliana e a Data Juliana Modificada
% Feito para a atividade de aula de GNC (01/07/2026)

    % Se não passar as horas/minutos/segundos, assume zero
    if nargin < 4
        hora = 0;
    end
    if nargin < 5
        minuto = 0;
    end
    if nargin < 6
        segundo = 0;
    end

    % variáveis auxiliares (I, J, K do algoritmo)
    I = ano;
    J = mes;
    K = dia;

    % Usando a função fix() para truncar a divisão inteira (igual ao FORTRAN)
    
    % dj1
    dj1 = K - 32075 + fix((1461 * (I + 4800 + fix((J - 14)/12))) / 4);
    
    % dj2
    dj2 = fix((367 * (J - 2 - 12 * fix((J - 14)/12))) / 12);
    
    % dj3
    dj3 = fix((-3 * fix((I + 4900 + fix((J - 14)/12)) / 100)) / 4);

    % data juliana final
    JD = dj1 + dj2 + dj3;

    % DJM0 (data juliana modificada sem fração, comeca a 0h UT)
    DJM0 = JD - 2400001;

    % calcula a fração do dia com base nas horas, minutos e segundos
    frac_dia = (hora + minuto/60 + segundo/3600) / 24;

    % DJM com a fração do dia
    DJM_frac = DJM0 + frac_dia;
end
