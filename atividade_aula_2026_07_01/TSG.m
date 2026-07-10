function tsg = TSG(ano, mes, dia, hora, minuto, segundo)
% Calcula o Tempo Sideral de Greenwich (theta_g) em graus
% Seção 5 do Escobal

    % DJ0: JD às 0h UT (-0.5 porque DiaJuliano referencia ao meio-dia)
    DJ0 = DiaJuliano(ano, mes, dia, 0, 0, 0) - 0.5;
    delta_t = hora * 3600 + minuto * 60 + segundo; % segundos do dia

    % Séculos Julianos desde Jan 0, 1900
    Tu = (DJ0 - 2415020.0) / 36525;

    % Dias no ano tropical
    L = 365.24219879 - 6.14e-6 * Tu;

    % Taxa de variação do TSG (graus/segundo)
    d_theta_dt = (1/240) * (1 + 1/L);

    % TSG às 0h UT
    theta_g0 = 99.69098329 + 36000.76893 * Tu + 3.87080e-4 * Tu^2;

    % TSG no instante pedido, reduzido para [0, 360)
    tsg = mod(theta_g0 + d_theta_dt * delta_t, 360);
end
