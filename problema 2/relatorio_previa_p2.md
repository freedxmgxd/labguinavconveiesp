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
   $$\epsilon = \frac{v^2}{2} - \frac{\mu}{r}$$
   $$a = -\frac{\mu}{2\epsilon}$$
   No MATLAB:
   ```matlab
   epsilon = (norm_v^2 / 2) - (mi / norm_r);
   a       = -mi / (2 * epsilon);
   ```
2. **Vetor excentricidade ($\vec{e}$) e magnitude ($e$):** Determina o achatamento da elipse e a direção do perigeu:
   $$\vec{e} = \frac{1}{\mu} \left[ \left(v^2 - \frac{\mu}{r}\right)\vec{r} - (\vec{r} \cdot \vec{v})\vec{v} \right]$$
   $$e = \|\vec{e}\|$$
   No MATLAB:
   ```matlab
   exc       = (1/mi) * ((norm_v^2 - mi/norm_r) * r - dot(r, v) * v);
   norma_exc = norm(exc);
   ```
3. **Momento angular específico ($\vec{h}$):** Perpendicular ao plano orbital:
   $$\vec{h} = \vec{r} \times \vec{v}$$
   $$h = \|\vec{h}\|$$
   No MATLAB:
   ```matlab
   h     = cross(r, v);
   h_mag = norm(h);
   ```
4. **Inclinação ($i$):** Ângulo de inclinação do plano da órbita em relação ao plano equatorial inercial ($\hat{K} = [0, 0, 1]^T$):
   $$i = \arccos\left(\frac{\vec{h} \cdot \hat{K}}{h}\right)$$
   No MATLAB:
   ```matlab
   inc   = acosd(dot([0 0 1], h) / h_mag);
   ```
5. **Linha dos nós e RAAN ($\Omega$):** Interseção do plano orbital com o plano equatorial:
   $$\vec{n} = \hat{K} \times \vec{h}$$
   $$\Omega = \arccos\left(\frac{\hat{I} \cdot \vec{n}}{\|\vec{n}\|}\right)$$
   onde $\hat{I} = [1, 0, 0]^T$. Se $n_y < 0$, então $\Omega = 360^\circ - \Omega$.
   
   No MATLAB:
   ```matlab
   n     = cross([0 0 1], h);
   ARNA  = acosd(dot([1 0 0], n) / norm(n));
   if n(2) < 0, ARNA = 360.0 - ARNA; end
   ```
   *A condicional corrige a ambiguidade da função arco-cosseno no plano equatorial.*
6. **Argumento do Perigeu ($\omega$):** Ângulo entre o nó ascendente e a direção do perigeu:
   $$\omega = \arccos\left(\frac{\vec{n} \cdot \vec{e}}{\|\vec{n}\| e}\right)$$
   onde se $e_z < 0$, então $\omega = 360^\circ - \omega$.
   
   No MATLAB:
   ```matlab
   OMEGA = acosd(dot(n, exc) / (norm(n) * norma_exc));
   if exc(3) < 0, OMEGA = 360.0 - OMEGA; end
   ```
7. **Anomalia Verdadeira inicial ($\nu_0$):** Posição angular do satélite em relação ao perigeu no instante inicial:
   $$\nu_0 = \arccos\left(\frac{\vec{e} \cdot \vec{r}}{e r}\right)$$
   onde se $\vec{r} \cdot \vec{v} < 0$, então $\nu_0 = 360^\circ - \nu_0$.
   
   No MATLAB:
   ```matlab
   NI = acosd(dot(exc, r) / (norma_exc * norm_r));
   if dot(r, v) < 0, NI = 360.0 - NI; end
   ```

8. **Parâmetros Derivados ($p, r_p, r_a, n, T$):**
   Semi-latus rectum, raio do perigeu, raio do apogeu, movimento médio e período orbital:
   $$p = \frac{h^2}{\mu}$$
   $$r_p = \frac{p}{1+e}$$
   $$r_a = \frac{p}{1-e}$$
   $$n = \sqrt{\frac{\mu}{a^3}}$$
   $$T = \frac{2\pi}{n}$$
   No MATLAB (no script `pb2.m`):
   ```matlab
   p     = h_mag^2 / mi;
   rp    = p / (1 + norma_exc);
   ra    = p / (1 - norma_exc);
   n_medio = sqrt(mi / a^3);
   T     = 2 * pi / n_medio;
   ```

### 1.3. Geração da Órbita no Plano Perifocal (Seção 4)
Para propagar a órbita inteira, cria-se um vetor contínuo de anomalias verdadeiras ($\nu \in [0, 2\pi]$) e calcula-se a distância polar $r$ e as coordenadas no plano perifocal $\vec{r}_p$:
$$r = \frac{p}{1 + e \cos\nu}$$
$$\vec{r}_p = \begin{bmatrix} r \cos\nu \\ r \sin\nu \\ 0 \end{bmatrix}$$

No MATLAB:
```matlab
ni_vetor = linspace(0, 2*pi, 1000);
r_polar  = p ./ (1 + norma_exc * cos(ni_vetor));
Ppolar = [r_polar .* cos(ni_vetor); r_polar .* sin(ni_vetor); zeros(1, length(ni_vetor))];
```
Isso gera coordenadas bidimensionais $(x_p, y_p)$ no referencial perifocal, onde o eixo $x_p$ aponta diretamente para o perigeu da órbita.

### 1.4. Transformação para o SGI com Rotação Horária (Seção 5)
Para transformar o vetor do plano perifocal para o SGI, realizam-se três rotações consecutivas utilizando ângulos negativos, pois o sistema perifocal está rotacionado de $(\Omega, i, \omega)$ em relação ao inercial. A matriz de rotação $Q$ é calculada pelas rotações elementares horárias:
$$R_{z1}(-\omega) = \begin{bmatrix} \cos(-\omega) & \sin(-\omega) & 0 \\ -\sin(-\omega) & \cos(-\omega) & 0 \\ 0 & 0 & 1 \end{bmatrix}$$
$$R_{x2}(-i) = \begin{bmatrix} 1 & 0 & 0 \\ 0 & \cos(-i) & \sin(-i) \\ 0 & -\sin(-i) & \cos(-i) \end{bmatrix}$$
$$R_{z3}(-\Omega) = \begin{bmatrix} \cos(-\Omega) & \sin(-\Omega) & 0 \\ -\sin(-\Omega) & \cos(-\Omega) & 0 \\ 0 & 0 & 1 \end{bmatrix}$$
$$Q = R_{z3}(-\Omega) R_{x2}(-i) R_{z1}(-\omega)$$
$$\vec{r}_{SGI} = Q \vec{r}_p$$

No MATLAB:
```matlab
rotz1 = [cosd(-OMEGA)  sind(-OMEGA) 0; -sind(-OMEGA) cosd(-OMEGA) 0; 0 0 1]; % Em torno do eixo Z do perifocal por -ω
rotx2 = [1 0 0; 0 cosd(-inc)  sind(-inc); 0 -sind(-inc) cosd(-inc)];         % Em torno do novo eixo X por -i
rotz3 = [cosd(-ARNA)   sind(-ARNA)  0; -sind(-ARNA)  cosd(-ARNA)  0; 0 0 1]; % Em torno do novo eixo Z por -Ω
Q     = rotz3 * rotx2 * rotz1;
Pxyz  = Q * Ppolar;
```
Adicionalmente, o código calcula o **Erro de Fechamento**, que mede a distância euclidiana entre a posição inicial fornecida pelo enunciado ($\vec{r}_{\text{inicial}}$) e o ponto calculado na órbita correspondente ao ângulo $\nu_0$ determinado. Para provar a exatidão matemática do algoritmo, calcula-se a posição exata analiticamente:
$$e_{\text{fechamento}} = \|\vec{r}_{\text{analítico}}(\nu_0) - \vec{r}_{\text{inicial}}\|_2$$

No MATLAB:
```matlab
Ppolar_exact = [p / (1 + norma_exc * cosd(NI)) * cosd(NI); p / (1 + norma_exc * cosd(NI)) * sind(NI); 0];
r_exact_sgi = Q * Ppolar_exact;
erro_fechamento = norm(r_exact_sgi' - r);
```

### 1.5. Latitude, Longitude e Efeito da Rotação Terrestre (Seção 6)
A projeção sub-satélite necessita da modelagem do movimento angular da Terra. As coordenadas esféricas (latitude e longitude) são obtidas a partir das coordenadas inerciais.

A latitude geocêntrica $\phi$ é dada por:
$$\phi = \arctan\left(\frac{r_z}{\sqrt{r_x^2 + r_y^2}}\right)$$

No MATLAB:
```matlab
lat = atan2d(rz, sqrt(rx.^2 + ry.^2));
```

Para a longitude, é necessário considerar o movimento de rotação da Terra. A partir do ângulo da anomalia verdadeira $\nu$, calcula-se o tempo decorrido desde a passagem pelo perigeu ($t$) usando as equações de Kepler. Com a formulação robusta e contínua em `atan2`, temos:
1. **Anomalia Excêntrica ($E$):** 
   $$E = 2 \arctan_2\left(\sqrt{1-e}\sin\left(\frac{\nu}{2}\right), \sqrt{1+e}\cos\left(\frac{\nu}{2}\right)\right)$$
2. **Anomalia Média ($M$) e Tempo ($t$):**
   $$M = E - e\sin(E) \quad \Rightarrow \quad t = \frac{M}{n}$$

No MATLAB:
```matlab
E_vetor = 2 * atan2(sqrt(1 - norma_exc) * sin(ni_vetor/2), sqrt(1 + norma_exc) * cos(ni_vetor/2));
M_vetor = E_vetor - norma_exc * sin(E_vetor);
t       = M_vetor / n_medio;
```

A longitude inercial é calculada e depois corrigida pelo efeito da rotação da Terra ($\omega_{\oplus}$) desde o instante inicial da simulação $t_{\text{sim}} = 0$ (instante em que o satélite passa pelo ponto de anomalia verdadeira inicial $\nu_0$, de forma que os eixos SGI e SCT coincidem perfeitamente nesse momento):
$$t_{\text{sim}} = t - t(\nu_0)$$
$$\lambda_{\text{inercial}} = \text{atan2}(r_y, r_x)$$
$$\lambda_E = \left( \lambda_{\text{inercial}} - \omega_{\oplus} t_{\text{sim}} + 180^\circ \pmod{360^\circ} \right) - 180^\circ$$

No MATLAB:
```matlab
t_sim = t - t(idx_ni);
lon = atan2d(ry, rx) - rad2deg(omega_terra * t_sim);
lon = mod(lon + 180, 360) - 180;
```
O código também substitui saltos instantâneos de longitude maiores que $180^\circ$ por `NaN` para evitar linhas espúrias cruzando o gráfico do mapa-múndi.

---

## 2. Análise dos Resultados

### 2.1. Elementos Orbitais Calculados
Ao executar o programa com os dados do enunciado, os seguintes resultados são exibidos no terminal (Command Window):

*   **Semieixo maior ($a$):** **$6732.6302 \text{ km}$** (indica uma órbita com altitude média próxima a $350\text{ km}$, típica de satélites em órbita baixa - LEO).
*   **Excentricidade ($e$):** **$0.001399$** (a órbita é praticamente circular, com desvio elíptico muito baixo).
*   **Inclinação ($i$):** **$42.7885^\circ$**.
*   **RAAN ($\Omega$):** **$159.6452^\circ$**.
*   **Arg. do Perigeu ($\omega$):** **$62.9484^\circ$**.
*   **Anomalia Verdadeira inicial ($\nu_0$):** **$38.9911^\circ$**.
*   **Período orbital ($T$):** **$5497.79 \text{ s}$** ($\approx 91.63 \text{ minutos}$).
*   **Erro de Fechamento (Analítico):** **$4.1926 \times 10^{-12} \text{ km}$** (erro nulo dentro da precisão de máquina, o que valida perfeitamente a consistência matemática das matrizes de transformação de coordenadas $Q$ na rotação 3D).

---

### 2.2. Discussão dos Gráficos Gerados

*   **Figura 1 (Plano da Órbita 2D - Cartesiano e Polar):** Exibe a órbita no plano bidimensional perifocal sob duas representações: uma cartesiana (sub-plot esquerdo) e outra polar real (sub-plot direito) em função da anomalia verdadeira $\nu$. Por conta da excentricidade baixíssima ($e \approx 0.0014$), ela se assemelha a uma circunferência perfeita centralizada na Terra (foco).
*   **Figura 2 (SGI 3D):** Mostra a orientação tridimensional da órbita em relação ao equador inercial. O plano da trajetória está inclinado a exatamente $42.79^\circ$. O vetor de posição inercial de entrada $\vec{R}$ está posicionado perfeitamente sobre a trajetória elíptica plotada.
*   **Figura 3 (Visualização da Geometria do Apêndice B.2):** Cria uma representação da esfera terrestre e ilustra a projeção da posição inercial $\vec{P}_{int}$ no plano equatorial, demarcando a latitude geocêntrica ($\Phi'$) e a longitude inercial ($\lambda_E$), funcionando como uma ilustração tridimensional fidedigna para a geometria descrita no referencial teórico.
*   **Figura 4 (Ground Track / Rastreio no Solo):** Plota o rastro da órbita sobre o mapa terrestre.
    1.  **Limites de Latitude:** A trajetória atinge uma latitude máxima de $42.79^\circ \text{ Norte}$ e mínima de $42.79^\circ \text{ Sul}$, condizente exatamente com a inclinação orbital real ($i = 42.79^\circ$).
    2.  **Deriva da Longitude:** O movimento no mapa não fecha em um ciclo fechado devido à rotação da Terra sob o veículo. O deslocamento (deriva para oeste) a cada órbita é de aproximadamente **$22.97^\circ$** (calculado via $\Delta\lambda = \omega_{\oplus} \cdot T = 7.292115\times10^{-5} \text{ rad/s} \times 5497.79 \text{ s} \approx 0.4009 \text{ rad} \approx 22.97^\circ$). O ponto inicial é plotado no local exato do rastreio correspondente a $t_{\text{sim}} = 0$.
