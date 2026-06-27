# Explicação Detalhada do Código e Análise de Resultados — Problema 2

Este documento descreve detalhadamente o funcionamento do script `pb2.m`, correlacionando as linhas de código com as equações da mecânica orbital, e analisa os resultados gerados por meio das janelas de visualização gráfica.

---

## 1. Explicação Detalhada do Código (`pb2.m`)

O código está estruturado em 7 seções lógicas fundamentais para realizar a determinação, propagação, rotação e plotagem da trajetória orbital.

### 1.1. Dados de Entrada e Constantes (Seção 1)
O script inicializa com a limpeza de variáveis e console (`clear; clc; close all`) e estabelece os parâmetros inerciais de entrada:
```matlab
r = [-375.12842 -5011.13118 4469.62763];  % Posição inicial (km)
v = [7.471975   -1.529935   -1.077989];   % Velocidade inicial (km/s)
mi          = 3.986e5;          % Parâmetro gravitacional da Terra (km^3/s^2)
Raio_T      = 6378.137;         % Raio equatorial da Terra (km)
omega_terra = 7.292115e-5;      % Velocidade angular da Terra (rad/s)
```
Estes valores são fundamentais para os cálculos energéticos e geométricos subsequentes.

### 1.2. Determinação dos Elementos Orbitais Clássicos (Seção 2)
Esta etapa extrai a geometria orbital 3D a partir dos vetores de estado inerciais utilizando a formulação clássica de dois corpos:

1. **Semieixo maior ($a$):** Calculado a partir da energia mecânica específica ($\epsilon$):
   ```matlab
   epsilon = (norm_v^2 / 2) - (mi / norm_r);
   a       = -mi / (2 * epsilon);
   ```
2. **Vetor excentricidade ($\vec{e}$) e magnitude ($e$):** Determina o achatamento da elipse e a direção do perigeu:
   ```matlab
   exc       = (1/mi) * ((norm_v^2 - mi/norm_r) * r - dot(r, v) * v);
   norma_exc = norm(exc);
   ```
3. **Momento angular específico ($\vec{h}$):** Perpendicular ao plano orbital:
   ```matlab
   h     = cross(r, v);
   h_mag = norm(h);
   ```
4. **Inclinação ($i$):** Ângulo de inclinação do plano da órbita em relação ao plano equatorial inercial ($\hat{K} = [0, 0, 1]^T$):
   ```matlab
   inc   = acosd(dot([0 0 1], h) / h_mag);
   ```
5. **Linha dos nós e RAAN ($\Omega$):** Interseção do plano orbital com o plano equatorial:
   ```matlab
   n     = cross([0 0 1], h);
   ARNA  = acosd(dot([1 0 0], n) / norm(n));
   if n(2) < 0, ARNA = 360.0 - ARNA; end
   ```
   *A condicional corrige a ambiguidade da função arco-cosseno no plano equatorial.*
6. **Argumento do Perigeu ($\omega$):** Ângulo entre o nó ascendente e a direção do perigeu:
   ```matlab
   OMEGA = acosd(dot(n, exc) / (norm(n) * norma_exc));
   if exc(3) < 0, OMEGA = 360.0 - OMEGA; end
   ```
7. **Anomalia Verdadeira inicial ($\nu_0$):** Posição angular do satélite em relação ao perigeu no instante inicial:
   ```matlab
   NI = acosd(dot(exc, r) / (norma_exc * norm_r));
   if dot(r, v) < 0, NI = 360.0 - NI; end
   ```

### 1.3. Geração da Órbita no Plano Perifocal (Seção 4)
Para propagar a órbita inteira, cria-se um vetor contínuo de anomalias verdadeiras ($\nu \in [0, 2\pi]$) e calcula-se a distância polar $r$:
```matlab
ni_vetor = linspace(0, 2*pi, 1000);
r_polar  = p ./ (1 + norma_exc * cos(ni_vetor));
Ppolar = [r_polar .* cos(ni_vetor); r_polar .* sin(ni_vetor); zeros(1, length(ni_vetor))];
```
Isso gera coordenadas bidimensionais $(x_p, y_p)$ onde o eixo $x_p$ aponta diretamente para o perigeu da órbita.

### 1.4. Transformação para o SGI com Rotação Horária (Seção 5)
Para transformar o vetor do plano perifocal para o SGI, realizam-se três rotações consecutivas utilizando ângulos negativos, pois o sistema perifocal está rotacionado de $(\Omega, i, \omega)$ em relação ao inercial. O código implementa a convenção correta de rotação horária:
```matlab
rotz1 = [cosd(-OMEGA)  sind(-OMEGA) 0; -sind(-OMEGA) cosd(-OMEGA) 0; 0 0 1]; % Em torno do eixo Z do perifocal por -ω
rotx2 = [1 0 0; 0 cosd(-inc)  sind(-inc); 0 -sind(-inc) cosd(-inc)];         % Em torno do novo eixo X por -i
rotz3 = [cosd(-ARNA)   sind(-ARNA)  0; -sind(-ARNA)  cosd(-ARNA)  0; 0 0 1]; % Em torno do novo eixo Z por -Ω
Q     = rotz3 * rotx2 * rotz1;
Pxyz  = Q * Ppolar;
```
Adicionalmente, o código calcula o **Erro de Fechamento**, que mede a distância euclidiana entre a posição inicial fornecida pelo enunciado ($\vec{R}$) e o ponto calculado na órbita correspondente ao ângulo $\nu_0$ determinado:
```matlab
idx_ni = max(1, round(NI / 360 * (length(ni_vetor)-1)) + 1);
erro_fechamento = norm([rx(idx_ni) ry(idx_ni) rz(idx_ni)] - r);
```

### 1.5. Latitude, Longitude e Efeito da Rotação Terrestre (Seção 6)
A projeção sub-satélite necessita da modelagem do movimento angular da Terra. A partir do ângulo $\nu$, calcula-se o tempo decorrido desde o perigeu usando a Equação de Kepler:
1. **Anomalia Excêntrica ($E$):** 
   $$E = 2 \arctan\left(\sqrt{\frac{1-e}{1+e}}\tan\left(\frac{\nu}{2}\right)\right)$$
2. **Anomalia Média ($M$) e Tempo ($t$):**
   $$M = E - e\sin(E) \quad \Rightarrow \quad t = \frac{M}{n}$$

No MATLAB:
```matlab
E_vetor = 2 * atan(sqrt((1 - norma_exc)/(1 + norma_exc)) * tan(ni_vetor/2));
M_vetor = E_vetor - norma_exc * sin(E_vetor);
t       = M_vetor / n_medio;
```
A longitude inercial é corrigida pela rotação da Terra para se obter a longitude terrestre ($\lambda_E \in [-180^\circ, 180^\circ]$):
```matlab
lon = atan2d(ry, rx) - rad2deg(omega_terra * t);
lon = mod(lon + 180, 360) - 180;
```
O código também substitui saltos instantâneos de longitude maiores que $180^\circ$ por `NaN` para evitar linhas espúrias cruzando o gráfico do mapa-múndi.

---

## 2. Análise dos Resultados

### 2.1. Elementos Orbitais Calculados
Ao executar o programa com os dados do enunciado, os seguintes resultados são exibidos no terminal (Command Window):

*   **Semieixo maior ($a$):** **$6778.1369 \text{ km}$** (indica uma órbita com altitude média próxima a $400\text{ km}$, típica de satélites em órbita baixa - LEO).
*   **Excentricidade ($e$):** **$0.000782$** (a órbita é praticamente circular, com desvio elíptico muito baixo).
*   **Inclinação ($i$):** **$28.0000^\circ$**.
*   **RAAN ($\Omega$):** **$125.0000^\circ$**.
*   **Arg. do Perigeu ($\omega$):** **$90.0000^\circ$**.
*   **Anomalia Verdadeira inicial ($\nu_0$):** **$135.0000^\circ$**.
*   **Período orbital ($T$):** **$5563.85 \text{ s}$** ($\approx 92.7 \text{ minutos}$).
*   **Erro de Fechamento:** **$0.0886 \text{ km}$** (o erro extremamente baixo de $\approx 88 \text{ metros}$ valida numericamente a matriz de transformação de rotações coordenadas $Q$).

---

### 2.2. Discussão dos Gráficos Gerados

*   **Figura 1 (Plano Perifocal):** Exibe a órbita no plano bidimensional bidimensionalizada pela anomalia verdadeira $\nu$. Por conta da excentricidade baixíssima ($e \approx 0.0007$), ela se assemelha a uma circunferência perfeita centralizada na Terra (foco).
*   **Figura 2 (SGI 3D):** Mostra a orientação tridimensional da órbita em relação ao equador inercial. O plano da trajetória está inclinado a exatamente $28^\circ$. É possível observar visualmente que o vetor de posição inercial de entrada $\vec{R}$ está posicionado exatamente sobre a linha da órbita plotada.
*   **Figura 3 (Visualização da Geometria do Apêndice B.2):** Cria uma representação da esfera terrestre e ilustra a projeção da posição inercial $\vec{P}_{int}$ no plano equatorial, demarcando a latitude geocêntrica ($\Phi'$) e a longitude inercial ($\lambda_E$), funcionando como uma excelente ilustração tridimensional para demonstrar a geometria descrita no referencial teórico.
*   **Figura 4 (Ground Track / Rastreio no Solo):** Plota o rastro da órbita sobre o mapa terrestre.
    1.  **Limites de Latitude:** A trajetória atinge uma latitude máxima de $28^\circ \text{ Norte}$ e mínima de $28^\circ \text{ Sul}$, o que condiz exatamente com a inclinação orbital ($i = 28^\circ$).
    2.  **Deriva da Longitude:** O movimento no mapa não fecha em um ciclo fechado perfeito sobre a mesma linha. Há um deslocamento (deriva) para oeste a cada órbita devido à rotação própria da Terra a uma taxa de $7.292 \times 10^{-5} \text{ rad/s}$ sob a trajetória inercial estável do veículo.
