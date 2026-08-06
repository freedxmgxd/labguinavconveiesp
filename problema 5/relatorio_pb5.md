# Relatório — Problema 5

**Curso:** Laboratório de Guiagem, Navegação e Controle de Veículos Espaciais
**Aluno:** Shinmen
**Tema:** Atitude de um V/E — Representação, Simulação e Determinação
**Código:** `pb5.m` (Apêndice A) — executar com `>> pb5`

---

## Sumário

1. [Descrição do problema](#1-descrição-do-problema)
2. [Introdução teórica](#2-introdução-teórica)
3. [Parte 1 — Simulação da rotação e recuperação dos ângulos](#3-parte-1--simulação-da-rotação-e-recuperação-dos-ângulos)
4. [Parte 2 — Determinação da atitude com sensor de estrelas](#4-parte-2--determinação-da-atitude-com-sensor-de-estrelas)
5. [Respostas às questões formuladas](#5-respostas-às-questões-formuladas)
6. [Pesquisa: sistemas inerciais, sensores de estrelas e Wahba](#6-pesquisa-sistemas-inerciais-sensores-de-estrelas-e-wahba)
7. [Conclusões](#7-conclusões)
8. [Referências](#8-referências)

---

## 1. Descrição do problema

O trabalho é composto por duas partes complementares, resolvidas por um único programa (`pb5.m`):

**Parte 1 — Problema direto e inverso da atitude.** Simula-se a rotação de um V/E em relação ao SGI, aplicando três rotações de Euler conhecidas ($\theta_{roll} = 187{,}9^\circ$, $\theta_{pitch} = 341{,}5^\circ$, $\theta_{yaw} = 79{,}7^\circ$) na sequência **132** a dois versores de referência $\hat r_1$ e $\hat r_2$. Obtidos os versores correspondentes no corpo, $\hat b_1$ e $\hat b_2$, aplica-se o **método algébrico (TRIAD)** para recuperar a matriz de atitude $A_T$ e, dela, extrair de volta os três ângulos. O sucesso do exercício é a recuperação exata dos ângulos simulados.

**Parte 2 — Problema real de determinação.** Um sensor rastreador estelar identifica duas estrelas (Matar/η Pegasi e Markab/α Pegasi) em seu campo de visada. Dadas as coordenadas inerciais catalogadas (AR/Dec, J2000.0) e as direções medidas no sistema do corpo, determina-se a atitude do V/E no instante da aquisição, em matriz de atitude e em ângulos de rolamento, arfagem e guinada.

---

## 2. Introdução teórica

### 2.1. O que é atitude

A **atitude** de um veículo espacial é a sua *orientação angular* no espaço, isto é, a relação entre o sistema de eixos fixo no corpo do V/E ($x$ = rolamento, $y$ = arfagem, $z$ = guinada) e um sistema de referência adotado. Quando o sistema de referência é inercial (SGI/J2000.0), fala-se em **atitude inercial**.

Atitude é um problema de **três graus de liberdade rotacionais**, completamente separado dos três graus de liberdade translacionais (posição/órbita). Determinar a atitude é responder: *"para onde o veículo está apontando?"* — não *"onde ele está?"*.

### 2.2. Representações da atitude

| Representação | Parâmetros | Vantagens | Desvantagens |
|---|---|---|---|
| **Matriz de cossenos diretores (DCM)** | 9 (6 vínculos) | Sem singularidades; composição direta por produto matricial | Redundante; custo de armazenamento e de reortonormalização |
| **Ângulos de Euler** | 3 | Interpretação física imediata (rolamento, arfagem, guinada) | Singularidade (*gimbal lock*); depende da sequência; trigonometria cara |
| **Quatérnios** | 4 (1 vínculo) | Sem singularidade; propagação numérica eficiente e estável | Não intuitivos; ambiguidade de sinal ($q$ e $-q$ representam a mesma atitude) |
| **Eixo/ângulo de Euler** | 4 ($\hat e$, $\Phi$) | Interpretação geométrica direta | Indefinido para $\Phi = 0$ |
| **Parâmetros de Rodrigues (Gibbs) / MRP** | 3 | Mínimo número de parâmetros | Singular em $\Phi = 180^\circ$ (Gibbs) ou $360^\circ$ (MRP) |

#### Matriz de cossenos diretores

A matriz de atitude $A$ é a matriz cujos elementos são os cossenos dos ângulos entre os eixos dos dois sistemas:

$$A_{ij} = \hat b_i \cdot \hat r_j$$

Ela transforma um vetor do sistema de referência para o sistema do corpo: $\vec v_b = A\,\vec v_r$. É uma matriz **ortogonal própria**: $A A^T = I$ e $\det(A) = +1$. Isso implica que a rotação **preserva normas e ângulos** — propriedade usada adiante como teste de consistência dos dados.

#### Ângulos de Euler e a matriz $A_{132}$

Adotando a convenção do curso, as três rotações elementares são:

$$R_x(\theta_{roll}) = \begin{bmatrix} 1 & 0 & 0 \\ 0 & c_r & s_r \\ 0 & -s_r & c_r \end{bmatrix},\quad
R_y(\theta_{pitch}) = \begin{bmatrix} c_p & 0 & -s_p \\ 0 & 1 & 0 \\ s_p & 0 & c_p \end{bmatrix},\quad
R_z(\theta_{yaw}) = \begin{bmatrix} c_y & s_y & 0 \\ -s_y & c_y & 0 \\ 0 & 0 & 1 \end{bmatrix}$$

A designação **132** indica a *ordem de aplicação* das rotações: primeiro em torno do eixo 1 ($x$), depois do eixo 3 ($z$), e por fim do eixo 2 ($y$). Como a rotação aplicada primeiro fica à direita no produto matricial:

$$\boxed{A_{132} = R_y(\theta_{pitch})\, R_z(\theta_{yaw})\, R_x(\theta_{roll}) = M_y M_z M_x}$$

Desenvolvendo o produto simbolicamente:

$$A_{132} = \begin{bmatrix}
c_p c_y & c_p s_y c_r + s_p s_r & c_p s_y s_r - s_p c_r \\
-s_y & c_y c_r & c_y s_r \\
s_p c_y & s_p s_y c_r - c_p s_r & s_p s_y s_r + c_p c_r
\end{bmatrix}$$

**Extração dos ângulos.** O elemento $A_{21} = -\sin\theta_{yaw}$ está isolado, o que dá o primeiro ângulo diretamente. Com ele conhecido, os outros dois saem de `atan2`, que resolve o quadrante corretamente:

$$\theta_{yaw} = \arcsin(-A_{21}), \qquad
\theta_{pitch} = \text{atan2}\!\left(\frac{A_{31}}{c_y}, \frac{A_{11}}{c_y}\right), \qquad
\theta_{roll} = \text{atan2}\!\left(\frac{A_{23}}{c_y}, \frac{A_{22}}{c_y}\right)$$

Limita-se o primeiro ângulo a $[-90^\circ, +90^\circ]$ para fugir da ambiguidade do seno. A singularidade (*gimbal lock*) ocorre quando $\cos\theta_{yaw} = 0$, ou seja, $\theta_{yaw} = \pm 90^\circ$: nesse caso $A_{11} = A_{31} = A_{22} = A_{23} = 0$ e rolamento e arfagem deixam de ser separáveis — só a sua soma/diferença é observável. O código trata esse caso explicitamente.

#### Quatérnios

O quatérnio de atitude codifica o eixo de Euler $\hat e$ e o ângulo de rotação $\Phi$:

$$\vec q = \begin{bmatrix} \hat e \sin(\Phi/2) \\ \cos(\Phi/2)\end{bmatrix} = \begin{bmatrix} q_1 & q_2 & q_3 & q_4 \end{bmatrix}^T, \qquad \|\vec q\| = 1$$

A matriz de atitude correspondente é

$$A(\vec q) = (q_4^2 - \vec q_v^T \vec q_v) I + 2 \vec q_v \vec q_v^T - 2 q_4 [\vec q_v \times]$$

Os quatérnios são o padrão em software de bordo por não terem singularidades e por a propagação cinemática ser **linear** em $\vec q$. Aparecem neste trabalho como o parâmetro natural do método q de Davenport (Seção 6.3).

### 2.3. Determinação da atitude: TRIAD e Wahba/QUEST

#### Método algébrico — TRIAD (Wertz, 1978)

Dados dois vetores não paralelos medidos no corpo ($\hat b_1, \hat b_2$) e conhecidos na referência ($\hat r_1, \hat r_2$), constrói-se uma **tríade ortonormal** em cada sistema:

$$\hat v_1 = \hat r_1, \quad \hat v_2 = \frac{\hat r_1 \times \hat r_2}{|\hat r_1 \times \hat r_2|}, \quad \hat v_3 = \hat v_1 \times \hat v_2 \;\Rightarrow\; M_r = [\hat v_1 \;|\; \hat v_2 \;|\; \hat v_3]$$

$$\hat w_1 = \hat b_1, \quad \hat w_2 = \frac{\hat b_1 \times \hat b_2}{|\hat b_1 \times \hat b_2|}, \quad \hat w_3 = \hat w_1 \times \hat w_2 \;\Rightarrow\; M_b = [\hat w_1 \;|\; \hat w_2 \;|\; \hat w_3]$$

Como ambas as tríades descrevem a *mesma* geometria física vista de dois sistemas, e são ortogonais por construção ($M_r^{-1} = M_r^T$):

$$\boxed{A_T = M_b\, M_r^T}$$

**Característica importante:** o TRIAD é *assimétrico*. Ele casa $\hat b_1$ **exatamente** e usa de $\hat b_2$ apenas a direção do produto vetorial — todo o erro de medida se acumula no segundo vetor. Por isso o vetor mais preciso deve ser sempre colocado na posição 1.

#### Problema de Wahba e o método q de Davenport (QUEST)

Com $N > 2$ direções medidas, o problema passa a ser sobredeterminado e resolve-se o **problema de Wahba (1965)**: achar a matriz ortogonal própria $A$ que minimiza

$$J(A) = \frac{1}{2}\sum_{k=1}^{N} w_k \left\| \hat b_k - A\,\hat r_k \right\|^2$$

com pesos $w_k$ inversamente proporcionais à variância de cada medida. Davenport mostrou que, escrevendo $A$ em quatérnios, o problema vira a maximização da forma quadrática $\vec q^T K \vec q$, com

$$B = \sum_k w_k\, \hat b_k \hat r_k^T, \qquad S = B + B^T, \qquad \sigma = \text{tr}(B), \qquad
\vec z = \begin{bmatrix} B_{23} - B_{32} \\ B_{31} - B_{13} \\ B_{12} - B_{21}\end{bmatrix}$$

$$K = \begin{bmatrix} S - \sigma I & \vec z \\ \vec z^{\,T} & \sigma \end{bmatrix} \quad (4\times 4)$$

A solução é o **autovetor de $K$ associado ao maior autovalor**. O **QUEST** (Shuster & Oh, 1981) é a implementação rápida desse resultado, que evita a decomposição completa em autovalores.

---

## 3. Parte 1 — Simulação da rotação e recuperação dos ângulos

### 3.1. Item 4 — Versores de referência

| Vetor | Componentes | Norma | Versor |
|---|---|---|---|
| $\vec r_1$ | $(1,\,-2,\,-2)$ | $3{,}0000$ | $\hat r_1 = (0{,}333333,\; -0{,}666667,\; -0{,}666667)$ |
| $\vec r_2$ | $(-1,\,2,\,3)$ | $3{,}7417$ | $\hat r_2 = (-0{,}267261,\; 0{,}534522,\; 0{,}801784)$ |

Ângulo entre $\hat r_1$ e $\hat r_2$: **168,5095°** (vetores quase antiparalelos, mas linearmente independentes — o produto vetorial não se anula, e o TRIAD é aplicável).

### 3.2. Item 5 — Matrizes de rotação

Com $\theta_{roll} = 187{,}9^\circ$, $\theta_{pitch} = 341{,}5^\circ$, $\theta_{yaw} = 79{,}7^\circ$:

$$M_x = \begin{bmatrix} 1 & 0 & 0 \\ 0 & -0{,}990509 & -0{,}137445 \\ 0 & 0{,}137445 & -0{,}990509 \end{bmatrix}\;
M_y = \begin{bmatrix} 0{,}948324 & 0 & 0{,}317305 \\ 0 & 1 & 0 \\ -0{,}317305 & 0 & 0{,}948324 \end{bmatrix}\;
M_z = \begin{bmatrix} 0{,}178802 & 0{,}983885 & 0 \\ -0{,}983885 & 0{,}178802 & 0 \\ 0 & 0 & 1 \end{bmatrix}$$

### 3.3. Item 6 — Versores no sistema do corpo

$$A_{132}^{sim} = M_y M_z M_x = \begin{bmatrix}
0{,}169562 & -0{,}880575 & -0{,}442535 \\
-0{,}983885 & -0{,}177105 & -0{,}024575 \\
-0{,}056735 & 0{,}439570 & -0{,}896415
\end{bmatrix}$$

$$\hat b_1 = A\,\hat r_1 = (0{,}938594,\; -0{,}193508,\; 0{,}285651) \qquad
\hat b_2 = A\,\hat r_2 = (-0{,}870822,\; 0{,}148583,\; -0{,}468607)$$

O ângulo entre $\hat b_1$ e $\hat b_2$ resulta **168,5095°**, idêntico ao do sistema de referência — confirmando que a rotação preserva os ângulos.

### 3.4. Item 7 — Método TRIAD

$$M_r = \begin{bmatrix}
0{,}333333 & -0{,}894427 & -0{,}298142 \\
-0{,}666667 & -0{,}447214 & 0{,}596285 \\
-0{,}666667 & 0 & -0{,}745356
\end{bmatrix} \qquad
M_b = \begin{bmatrix}
0{,}938594 & 0{,}242144 & -0{,}245781 \\
-0{,}193508 & 0{,}959217 & 0{,}206050 \\
0{,}285651 & -0{,}145837 & 0{,}947172
\end{bmatrix}$$

$$\boxed{A_T = M_b M_r^T = \begin{bmatrix}
0{,}169562 & -0{,}880575 & -0{,}442535 \\
-0{,}983885 & -0{,}177105 & -0{,}024575 \\
-0{,}056735 & 0{,}439570 & -0{,}896415
\end{bmatrix}}$$

**Verificações:**

| Teste | Resultado |
|---|---|
| Ortogonalidade $\|A_T A_T^T - I\|$ | $4{,}45 \times 10^{-16}$ |
| $\det(A_T)$ | $1{,}0000000000$ |
| Erro vs. matriz simulada $\|A_T - A_{132}^{sim}\|$ | $3{,}35 \times 10^{-16}$ |

A matriz determinada pelo TRIAD é **idêntica** à matriz simulada, a menos do erro de arredondamento de ponto flutuante.

### 3.5. Item 8 — Ângulos recuperados

| Ângulo | Simulado | Recuperado | Equivalente em $[0^\circ,360^\circ)$ | Erro |
|---|---|---|---|---|
| Rolamento (*roll*) | $187{,}9000^\circ$ | $-172{,}1000^\circ$ | $187{,}9000^\circ$ | $0{,}0 \times 10^{0}$ |
| Arfagem (*pitch*) | $341{,}5000^\circ$ | $-18{,}5000^\circ$ | $341{,}5000^\circ$ | $5{,}7 \times 10^{-14}$ |
| Guinada (*yaw*) | $79{,}7000^\circ$ | $79{,}7000^\circ$ | $79{,}7000^\circ$ | $-2{,}8 \times 10^{-14}$ |

Os três ângulos foram **recuperados exatamente**. A função `atan2` devolve o valor principal em $(-180^\circ, 180^\circ]$, de modo que $187{,}9^\circ$ e $341{,}5^\circ$ reaparecem como $-172{,}1^\circ$ e $-18{,}5^\circ$: são o **mesmo ângulo** (diferença de exatamente $360^\circ$) e geram rigorosamente a mesma matriz de atitude. O fechamento do ciclo confirma isso:

$$\left\| A_{132}(\theta_{recuperados}) - A_T \right\| = 4{,}13 \times 10^{-16}$$

Isso responde à observação do enunciado: os três ângulos relacionam o sistema de referência com o sistema do corpo e, portanto, **representam a atitude do V/E em ângulos de Euler**.

---

## 4. Parte 2 — Determinação da atitude com sensor de estrelas

### 4.1. Conversão das coordenadas catalogadas

$$\alpha[^\circ] = \left(h + \frac{m}{60} + \frac{s}{3600}\right)\times 15 \qquad
\hat r = \begin{bmatrix} \cos\delta\cos\alpha \\ \cos\delta\sin\alpha \\ \sin\delta \end{bmatrix}$$

| Estrela | AR | Dec | $\alpha$ [°] | $\delta$ [°] | Versor no SGI |
|---|---|---|---|---|---|
| **Matar** (η Pegasi) | 22h 44m 16,0s | +30° 21′ 32,7″ | 341,066667 | 30,359083 | $(0{,}816190,\; -0{,}279975,\; 0{,}505418)$ |
| **Markab** (α Pegasi) | 23h 06m 06,0s | +15° 20′ 53,8″ | 346,525000 | 15,348278 | $(0{,}937788,\; -0{,}224710,\; 0{,}264686)$ |

### 4.2. Verificação de consistência dos dados — achado importante

Antes de determinar a atitude, aplica-se o teste obrigatório: **uma rotação preserva o ângulo entre vetores**, logo o ângulo entre $\vec b_1$ e $\vec b_2$ tem de ser igual ao ângulo entre $\hat r_1$ e $\hat r_2$.

| Grandeza | Valor |
|---|---|
| Separação angular Matar–Markab no SGI | **15,8240°** |
| Ângulo entre $\vec b_1$ e $\vec b_2$ conforme publicados | **54,9466°** |
| **Discrepância** | **+39,1227°** |

Um sensor rastreador estelar tem exatidão da ordem de **segundos de arco**. Uma discrepância de 39° não é ruído de medida: **os dados publicados são geometricamente inconsistentes** — não existe matriz de rotação alguma que leve $\hat r_1 \to \vec b_1$ e $\hat r_2 \to \vec b_2$ simultaneamente.

Testando a hipótese de erro de sinal na primeira componente de $\vec b_2$:

| $\vec b_2$ | Ângulo $\angle(\vec b_1, \vec b_2)$ | Discrepância |
|---|---|---|
| $[-0{,}3521\;\; -0{,}9351\;\; 0{,}0386]$ (publicado) | 54,9466° | +39,1227° |
| $[+0{,}3521\;\; -0{,}9351\;\; 0{,}0386]$ (sinal corrigido) | **15,8238°** | **−0,0002° (0,7″)** |

Com o sinal invertido, a concordância com o catálogo é de **0,7 segundos de arco** — exatamente o resíduo esperado do arredondamento das componentes em 4 casas decimais. Conclui-se que **o sinal “−” em $b_2(1)$ é um erro de digitação do enunciado**. O programa resolve e reporta os dois cenários.

### 4.3. Resultados

**Cenário A — dados exatamente como publicados**

$$A = \begin{bmatrix}
0{,}017882 & -0{,}123381 & 0{,}992198 \\
-0{,}999170 & 0{,}034116 & 0{,}022250 \\
-0{,}036595 & -0{,}991773 & -0{,}122669
\end{bmatrix}$$

**Cenário B — sinal de $b_2(1)$ corrigido (adotado)**

$$\boxed{A = \begin{bmatrix}
0{,}160587 & 0{,}232954 & 0{,}959137 \\
-0{,}952956 & 0{,}289687 & 0{,}089193 \\
-0{,}257072 & -0{,}928339 & 0{,}268515
\end{bmatrix}}$$

**Atitude em ângulos de Euler (sequência 132):**

| Ângulo | Cenário A | **Cenário B (adotado)** | Equivalente em $[0^\circ,360^\circ)$ |
|---|---|---|---|
| Rolamento (*roll*) | 33,1118° | **17,1133°** | 17,1133° |
| Arfagem (*pitch*) | −63,9579° | **−58,0079°** | 301,9921° |
| Guinada (*yaw*) | 87,6657° | **72,3557°** | 72,3557° |
| Resíduo angular em $\vec b_2$ | **39,1227°** | **0,0002°** | — |

Em ambos os cenários a matriz é ortogonal ($\|AA^T - I\| \sim 10^{-16}$, $\det A = +1$) e o resíduo em $\vec b_1$ é **exatamente zero** — consequência direta da assimetria do TRIAD. A diferença está no segundo vetor: no cenário A sobra um erro de 39°, denunciando a inconsistência; no cenário B o resíduo cai para 0,0002° (menos de 1″), compatível com a precisão dos dados.

**O cenário B é a resposta adotada**, por ser o único fisicamente admissível.

---

## 5. Respostas às questões formuladas

### Questão 1 — Atitude inercial: o que é? Quantas e quais maneiras há de representá-la?

**Atitude inercial** é a orientação angular do sistema de eixos fixo no corpo do V/E em relação a um **sistema de referência inercial** — aqui o SGI/J2000.0 (geocêntrico equatorial, com $X$ apontando para o equinócio vernal e $Z$ para o polo norte celeste). É a atitude "absoluta", medida contra um referencial que não gira, e é exatamente o que um sensor de estrelas mede: as estrelas são as marcas fixas desse referencial.

As representações usuais estão tabeladas na Seção 2.2. Resumidamente, há **cinco famílias principais**: matriz de cossenos diretores (9 parâmetros), ângulos de Euler (3), quatérnios (4), eixo/ângulo de Euler (4) e parâmetros de Rodrigues/MRP (3). Nenhuma representação de apenas 3 parâmetros é globalmente não singular — esse é um resultado topológico, e é a razão de os quatérnios dominarem o software de bordo.

**A matriz de cossenos diretores** é aquela cujos elementos são os cossenos dos ângulos entre os eixos dos dois sistemas, $A_{ij} = \hat b_i \cdot \hat r_j$; suas linhas são os versores do corpo escritos no sistema de referência. **As matrizes de atitude** são ortogonais próprias ($A A^T = I$, $\det A = +1$), formando o grupo $SO(3)$: preservam comprimentos, ângulos e orientação, e sua inversa é a transposta.

**A matriz $A_{132}$** é a matriz de atitude construída aplicando-se as rotações elementares na ordem $x \to z \to y$, ou seja, $A_{132} = M_y M_z M_x$ — sua forma simbólica e as fórmulas de extração dos ângulos estão na Seção 2.2, e ela é a matriz usada em ambas as partes deste trabalho.

**Contraste com a matriz $A_{231}$** (citada nas referências do curso): a sequência 231 aplica as rotações na ordem $y \to z \to x$, resultando em $A_{231} = R_x R_z R_y$, com forma simbólica

$$A_{231} = \begin{bmatrix}
c_p c_y & s_y & -c_y s_p \\
s_p s_r - c_p c_r s_y & c_r c_y & c_p s_r + c_r s_p s_y \\
c_r s_p + c_p s_r s_y & -c_y s_r & c_p c_r - s_p s_r s_y
\end{bmatrix}$$

de onde $\theta_{yaw} = \arcsin(A_{12})$, $\theta_{pitch} = \text{atan2}(A_{13}/-c_y,\, A_{11}/c_y)$ e $\theta_{roll} = \text{atan2}(A_{32}/-c_y,\, A_{22}/c_y)$. O padrão de extração é sempre o mesmo: **localizar na matriz simbólica o elemento que isola um único ângulo** (aqui $A_{12} = s_y$; na 132, $A_{21} = -s_y$), extraí-lo com o arco-seno limitado a $[-90^\circ, 90^\circ]$, e obter os outros dois com `atan2` usando o primeiro já conhecido.

### Questão 2 — Há problema na suposição do enunciado, já que a origem do SGI está no centro da Terra e a origem do sistema do corpo é o CG do veículo?

**Para a determinação da atitude, não há problema algum.** A atitude é uma propriedade **puramente rotacional**: depende apenas da orientação relativa entre as *direções* dos eixos, e uma translação da origem não altera a direção de eixo nenhum. Matematicamente, a matriz $A$ relaciona versores, e versores são invariantes por translação. Por isso é lícito e usual "transportar" o SGI para o CG do veículo — o sistema resultante é paralelo ao SGI e serve de referência para a atitude.

**Onde a diferença de origens efetivamente importa:**

1. **Paralaxe.** Vetores medidos para alvos *próximos* dependem da posição. Para **estrelas**, no entanto, a distância é tão grande que o efeito é nulo na prática: a estrela mais próxima (Proxima Centauri, 1,3 pc) tem paralaxe de 0,77″ para uma base de 1 UA; para uma base de ~7000 km (raio de uma órbita baixa) a paralaxe cai para ~$4\times10^{-5}$ segundos de arco — dez mil vezes menor que a resolução de qualquer sensor. **É justamente por isso que o sensor de estrelas é o instrumento de atitude por excelência.** Já para sensores de Terra, de Sol ou de Lua o efeito de paralaxe é significativo e precisa ser corrigido com o conhecimento da órbita.

2. **Aberração da luz.** Este é o efeito que *realmente* acopla órbita e atitude num sensor estelar de alta precisão. A velocidade do observador desloca a direção aparente da estrela em $\Delta\theta \approx v/c$:
   - aberração anual (Terra a 29,8 km/s): **~20,5″**;
   - aberração devida à velocidade orbital do satélite (LEO, ~7,5 km/s): **~5,2″**.

   Ambas são muito maiores que a exatidão de um AST (poucos segundos de arco) e devem ser compensadas usando a posição e a velocidade orbitais — ou seja, **atitude e órbita se acoplam pela aberração, não pela translação da origem**.

3. **Acoplamento dinâmico.** Na dinâmica (não na cinemática), a escolha do CG como origem do sistema do corpo é o que desacopla translação de rotação nas equações de Euler. Se a origem não fosse o CG, apareceriam termos cruzados de inércia. Além disso, o CG **migra** durante a missão (consumo de propelente, deslocamento de apêndices), e o desalinhamento resultante entre o eixo de empuxo e o CG produz torques de perturbação que o SCA precisa absorver.

**Conclusão:** a suposição do enunciado é válida para o propósito do problema. A ressalva técnica correta não é a translação da origem, e sim a **aberração estelar**, que num sensor real de arcossegundos é obrigatória.

### Questão 3 — Qual a diferença entre atitude (inercial) e desvios em atitude?

São conceitos de **natureza e magnitude diferentes**, frequentemente confundidos porque ambos são expressos em rolamento/arfagem/guinada:

| | **Atitude (inercial)** | **Desvios em atitude** |
|---|---|---|
| **Referência** | Sistema inercial (SGI/J2000.0) | Atitude *nominal/desejada* (p. ex. o sistema orbital LVLH, ou uma atitude comandada) |
| **Magnitude** | Ângulos arbitrários, $0$ a $360^\circ$ | Ângulos **pequenos**, tipicamente frações de grau |
| **Natureza** | Grandeza **absoluta** — a orientação em si | Grandeza **diferencial** — o *erro* de apontamento |
| **Matemática** | Exige a matriz completa; não linearizável | Permite linearização: $\sin\varepsilon \approx \varepsilon$, $\cos\varepsilon \approx 1$, e $A \approx I - [\vec\varepsilon\times]$ |
| **Obtida por** | **Determinação** de atitude (sensores + TRIAD/QUEST) | Diferença entre atitude determinada e atitude de referência |
| **Usada em** | Navegação, conhecimento absoluto de apontamento | **Controle** — é o sinal de erro que realimenta o SCA |

Formalmente, se $A_{det}$ é a atitude determinada e $A_{ref}$ a de referência, o desvio é a rotação residual

$$\delta A = A_{det}\, A_{ref}^{T} \approx I - [\vec\varepsilon \times]$$

e $\vec\varepsilon = (\varepsilon_{roll}, \varepsilon_{pitch}, \varepsilon_{yaw})$ é o vetor de desvios. É a pequenez de $\vec\varepsilon$ que permite **linearizar** a dinâmica e projetar controladores lineares (PID, LQR) para o sistema de controle de atitude — o que seria impossível com a atitude absoluta, cuja cinemática é essencialmente não linear.

Em resumo: **a atitude diz para onde o veículo aponta; os desvios dizem o quanto ele está errando em relação a para onde deveria apontar.** Um satélite de observação da Terra tem atitude inercial que varia continuamente ao longo da órbita (ele acompanha o nadir), enquanto seus desvios em atitude permanecem próximos de zero — é exatamente esse o trabalho do SCA.

---

## 6. Pesquisa: sistemas inerciais, sensores de estrelas e Wahba

### 6.1. Sistemas de coordenadas inerciais

Um referencial rigorosamente inercial não existe; na prática usam-se referenciais **quase-inerciais**, definidos por convenção e materializados por catálogos de objetos distantes.

| Sistema | Origem | Eixo $X$ | Plano fundamental | Uso |
|---|---|---|---|---|
| **GCRF / J2000.0 (SGI)** | Centro da Terra | Equinócio vernal médio em 01/01/2000, 12h TT | Equador médio de J2000.0 | **Padrão em dinâmica orbital e de atitude de satélites terrestres** |
| **ICRF** | Baricentro do Sistema Solar | Fixo por ~300 radiofontes extragalácticas (quasares) | Próximo ao equador J2000 | Definição fundamental moderna (IAU, 1998); base de todos os demais |
| **BCRF / ICRS baricêntrico** | Baricentro do Sistema Solar | Idem ICRF | — | Efemérides planetárias, missões interplanetárias |
| **MOD / TOD** | Centro da Terra | Equinócio médio / verdadeiro da data | Equador da data | Etapas intermediárias na redução precessão–nutação |
| **Eclíptico J2000** | Sol ou baricentro | Equinócio vernal J2000 | Eclíptica | Trajetórias interplanetárias |

**Diferença essencial entre J2000 e ICRF:** o J2000 é definido *dinamicamente* (pelo movimento da Terra — equador e equinócio de uma época), enquanto o ICRF é definido *cinematicamente* (por quasares, sem referência a qualquer movimento). O ICRF é hoje a realização fundamental; o GCRF é sua contraparte geocêntrica, e difere do J2000 dinâmico por apenas ~$0{,}02''$ (o *frame bias*), diferença desprezível para a maioria das aplicações mas relevante em geodésia de precisão.

**Sistema utilizado neste estudo:** as coordenadas das estrelas foram fornecidas em **AR/Dec no SGI equatorial geocêntrico, época J2000.0**, que é o sistema em que os catálogos estelares modernos (Hipparcos, Tycho-2, Gaia) publicam suas posições. Vale notar que catálogos de alta precisão fornecem também **movimento próprio**, e a redução rigorosa de posições catalogadas para direções aparentes na época da observação exige aplicar movimento próprio, paralaxe, precessão–nutação e aberração.

### 6.2. Sensores de estrelas

**O que são.** São câmeras digitais de alta sensibilidade que fotografam um trecho do céu, identificam as estrelas presentes na imagem comparando-as com um catálogo embarcado, e daí calculam a atitude inercial do V/E. São os sensores de atitude **mais precisos** disponíveis, e os únicos que fornecem atitude inercial *absoluta e de três eixos* a partir de uma única medida.

**Como funcionam** — cadeia de processamento típica:

1. **Aquisição** — a óptica projeta o campo de visada (CDV, tipicamente 8°×8° a 20°×20°) sobre um detector CCD ou APS/CMOS. O foco é propositalmente levemente *desfocado*, para espalhar a estrela por vários pixels.
2. **Detecção e centroidagem** — identificam-se as manchas luminosas acima do limiar de ruído e calcula-se o **centroide** (linha, coluna) de cada uma com precisão *subpixel* (~1/10 de pixel), justamente graças ao desfoque proposital. É esta a etapa que o enunciado da Parte 2 menciona.
3. **Conversão para direções** — com o modelo geométrico da câmera (distância focal, ponto principal, distorções), cada centroide vira um versor $\vec b_k$ no sistema fixo no corpo. **São exatamente os $\vec b_1$ e $\vec b_2$ do problema.**
4. **Identificação (*star ID*)** — o padrão geométrico de distâncias angulares entre as estrelas detectadas é comparado com o catálogo embarcado. No modo *lost-in-space* (sem atitude prévia) usam-se algoritmos de triângulos, *grids* ou *k-vector*; no modo de rastreio, a atitude anterior restringe drasticamente a busca.
5. **Determinação da atitude** — de posse dos pares $(\hat r_k, \vec b_k)$, aplica-se TRIAD (2 estrelas) ou QUEST/Wahba ($N$ estrelas), obtendo a matriz de atitude ou o quatérnio.

**Precisão e estado da arte.** Unidades modernas atingem **1–3″ (1σ)** nos eixos transversais e ~10–20″ em torno do eixo de visada (*boresight*) — a assimetria vem de a rotação em torno do eixo óptico ser observada por uma alavanca angular pequena (metade do CDV). Taxas de atualização de 4–10 Hz são comuns, com massa de 0,5–3 kg e consumo de 5–15 W. A tendência atual são sensores **autônomos** (com todo o processamento e o catálogo embarcados, entregando quatérnio pronto), miniaturizados para CubeSats, e com múltiplas cabeças ópticas para robustez.

**Limitações operacionais:** cegueira por luz espúria (Sol, Terra, Lua no CDV — daí os *baffles*), degradação por radiação e por *hot pixels*, e erros dinâmicos de arrasto de imagem sob altas taxas angulares.

**Observação do enunciado sobre estrelas brilhantes.** O texto assinala corretamente que estrelas muito brilhantes — como Sirius, Vega, Aldebaran, e também Matar e Markab — em geral **não** são usadas pelos sensores. Os motivos são: **saturação** dos pixels centrais, que destrói a informação de forma da mancha e portanto a centroidagem subpixel; *blooming* (vazamento de carga para pixels vizinhos em CCDs); e o fato de as estrelas mais brilhantes serem poucas e mal distribuídas no céu, gerando cobertura irregular. Catálogos de sensores estelares tipicamente trabalham na faixa de **magnitude 2 a 6**, onde há milhares de estrelas bem distribuídas e não saturantes.

**O AST/ADAST (INPE–OMNISYS–UFABC).** O *Autonomous Star Tracker* é o sensor de estrelas autônomo desenvolvido no âmbito da parceria INPE–OMNISYS com participação da UFABC, integrante da Plataforma Multi-Missão (PMM) brasileira. O **ADAST** é o software simulador de operação do instrumento, usado em desenvolvimento e testes. Trata-se de um avanço tecnológico nacional relevante, por dar autonomia ao país num componente crítico e de exportação controlada (de Brum et al., 2013).

### 6.3. Questão C — Uso de muitas estrelas no CDV

**Pergunta:** o sensor identifica frequentemente 8, 10, 12, 16 estrelas no CDV, mas o método algébrico usa apenas duas direções. Como se procede?

**Resposta:** não se usa o método algébrico. Passa-se do problema *determinado* (2 vetores, 3 incógnitas, solução exata) para um problema **sobredeterminado** ($N$ vetores), resolvido por **mínimos quadrados** — que é precisamente o **problema de Wahba** (Seção 2.3). Na prática:

1. **Descartar informação é inaceitável.** Usar só 2 de $N$ estrelas joga fora medidas boas e mantém o erro de uma única estrela dominando o resultado. Como os erros são aproximadamente independentes, a média ponderada de $N$ medidas reduz o erro por um fator $\approx 1/\sqrt{N}$.
2. **Ponderação.** Cada estrela recebe peso $w_k \propto 1/\sigma_k^2$: estrelas mais fracas ou próximas da borda do CDV (onde a distorção óptica é maior) contribuem menos.
3. **Solução ótima.** Monta-se a matriz $K$ de Davenport ($4\times4$) e toma-se o autovetor associado ao maior autovalor — o quatérnio de atitude. O **QUEST** obtém a mesma solução sem decomposição completa em autovalores, resolvendo a equação característica por Newton–Raphson a partir de $\lambda_0 \approx \sum w_k$, o que o torna adequado ao processamento de bordo em tempo real. Outras soluções: **SVD** (a mais robusta numericamente), **FOAM**, **ESOQ/ESOQ-2**.
4. **Vantagem adicional: detecção de erros.** Com redundância é possível testar a consistência de cada medida e **rejeitar** estrelas mal identificadas (*outliers*), o que é impossível com apenas duas direções.

**Demonstração numérica** (implementada em `pb5.m`, Seção C.2). Com uma atitude verdadeira conhecida, 8 estrelas sintéticas no CDV e ruído gaussiano de 20″ por eixo:

| Método | Direções usadas | Erro de atitude |
|---|---|---|
| **TRIAD** | 2 das 8 estrelas | **39,56″** |
| **QUEST / Wahba** | todas as 8 estrelas | **16,54″** |

O ganho de aproximadamente $\sqrt{N/2} = 2$ confirma a previsão teórica. Aplicado às duas estrelas reais do problema (cenário B), o método q de Davenport reproduz o resultado do TRIAD com diferença de apenas **0,0003°** — como esperado, pois com apenas duas direções quase consistentes ambos os métodos convergem para a mesma solução; a diferença só se manifesta com redundância e ruído.

---

## 7. Conclusões

1. **A Parte 1 fechou o ciclo completo com precisão de máquina.** Simuladas as três rotações na sequência 132 com $\theta_{roll} = 187{,}9^\circ$, $\theta_{pitch} = 341{,}5^\circ$ e $\theta_{yaw} = 79{,}7^\circ$, o método TRIAD recuperou a matriz de atitude com erro de $3{,}3\times10^{-16}$ e os três ângulos foram recuperados exatamente. Isso valida simultaneamente a implementação das matrizes de rotação, do método algébrico e das fórmulas de extração de ângulos da matriz $A_{132}$.

2. **A ambiguidade de $360^\circ$ é aparente, não real.** Os ângulos retornam como $-172{,}1^\circ$ e $-18{,}5^\circ$ porque `atan2` devolve o valor principal; são os mesmos ângulos de $187{,}9^\circ$ e $341{,}5^\circ$ e produzem rigorosamente a mesma matriz de atitude, como o teste de fechamento confirmou.

3. **A Parte 2 revelou uma inconsistência nos dados do enunciado.** O teste de invariância angular — que deve preceder toda determinação de atitude — mostrou que as direções medidas formam 54,95° entre si, contra 15,82° de separação real entre Matar e Markab no catálogo. Nenhuma rotação pode conciliar as duas medidas. Invertendo o sinal da primeira componente de $\vec b_2$, a concordância passa a ser de **0,7 segundos de arco**, identificando um erro de digitação no enunciado. **Atitude adotada: rolamento 17,1133°, arfagem −58,0079° (≡ 301,9921°), guinada 72,3557°**, com resíduo de 0,0002° na segunda estrela.

4. **O teste de consistência é parte do método, não um extra.** Este problema ilustra na prática por que sistemas reais de determinação de atitude sempre verificam a preservação dos ângulos antes de aceitar uma medida: o TRIAD produz *silenciosamente* uma matriz perfeitamente ortogonal, com $\det = +1$, mesmo a partir de dados fisicamente impossíveis. A matriz "parece" válida; só o resíduo de 39° denuncia o problema.

5. **A assimetria do TRIAD é sua principal limitação.** O método casa o primeiro vetor exatamente e despeja todo o erro no segundo. Daí a regra prática de colocar sempre a medida mais confiável na posição 1 — e a superioridade do QUEST/Wahba, que distribui o erro de forma ótima entre todas as observações, com ganho de $\approx 1/\sqrt{N}$ confirmado numericamente (39,56″ → 16,54″ com 8 estrelas).

6. **A diferença de origens entre SGI e sistema do corpo é irrelevante para a atitude**, por ser esta uma grandeza puramente rotacional, e porque a paralaxe estelar a partir de uma órbita terrestre é ~$10^{-5}$ segundos de arco. O efeito que realmente acopla órbita e atitude num sensor estelar de precisão é a **aberração da luz** (~5″ em órbita baixa, ~20″ de aberração anual), que excede a exatidão do instrumento e deve ser corrigida com o conhecimento da velocidade orbital.

---

## 8. Referências

- **WERTZ, J. R.** *Spacecraft Attitude Determination and Control.* London: D. Reidel, 1978.
- **SHUSTER, M. D.; OH, S. D.** "Three-axis attitude determination from vector observations". *Journal of Guidance and Control*, vol. 4, n. 1, 1981, pp. 70–77.
- **WAHBA, G.** "A Least Squares Estimate of Satellite Attitude". *SIAM Review*, vol. 7, n. 3, 1965, p. 409.
- **SIDI, M. J.** *Spacecraft Dynamics and Control: A Practical Engineering Approach.* Cambridge University Press, 1997.
- **de BRUM, A. G. V.; FIALHO, M. A.; SELINGARDI, M. L.; BORREGO, N. D.; LOURENÇO, J. L.** "The Brazilian Autonomous Star Tracker – AST". *WSEAS Transactions on Systems*, vol. 12, 2013, pp. 459–470. Disponível em: http://www.wseas.us/journal/pdf/systems/2013/a045702-262.pdf
- **de BRUM, A. G. V.** Notas de aula — ESTS003-17 Introdução à Astronáutica e Lab. de GN&C, CECS/UFABC.
- **de BRUM, A. G. V.** "Sobre outras matrizes de Atitude (312, 132, etc.)" — material de apoio da disciplina.
- **MARKLEY, F. L.; CRASSIDIS, J. L.** *Fundamentals of Spacecraft Attitude Determination and Control.* Springer, 2014.
- Coordenadas estelares: Wikipedia e software **Stellarium** (https://stellarium.org).

---

## Apêndice A — Código-fonte

O programa `pb5.m` resolve as duas partes e a questão C. Executar com `>> pb5` (compatível com MATLAB e GNU Octave). Saídas gráficas: `pb5_parte1_vetores.png` e `pb5_parte2_estrelas.png`.

Estrutura:

| Seção | Conteúdo |
|---|---|
| Parte 1, itens 4–5 | Versores de referência; matrizes $M_x$, $M_y$, $M_z$ |
| Parte 1, item 6 | $A_{132} = M_yM_zM_x$; geração de $\hat b_1$, $\hat b_2$ |
| Parte 1, itens 7–8 | TRIAD ($M_r$, $M_b$, $A_T$); extração e verificação dos ângulos |
| Parte 2 | Conversão AR/Dec → versores; **teste de consistência**; TRIAD nos dois cenários |
| Questão C | Método q de Davenport; simulação TRIAD × QUEST com 8 estrelas |
| Funções auxiliares | `rotx_gnc`, `roty_gnc`, `rotz_gnc`, `triad`, `ang_from_A132`, `davenport_q`, `hms2deg`, `dms2deg`, `radec2vec`, `err_ang_matriz` |
