# Relatório da Atividade - Tempo e Data Juliana
**Disciplina**: Laboratório de Guiagem, Navegação e Controle de Veículos Espaciais (UFABC)  
**Atividade**: atividade_aula_2026_07_01  

---

## 1. O que é Tempo Universal (UT) e Data Juliana (DJ)?

* **Tempo Universal (UT)**: É a medida do tempo baseada na rotação da Terra em relação a estrelas e outros corpos celestes. É a referência principal de tempo para a navegação e atividades espaciais.
* **Data Juliana (DJ)**: É um sistema que conta de forma contínua a quantidade de dias passados desde um instante de referência. Não usa meses ou anos, sendo apenas um número decimal/inteiro contínuo, facilitando os cálculos de órbitas e tempo espacial.

---

## 2. Perguntas do PDF

### A. Que data de referência é essa? O que ela tem de especial?
A referência é o **meio-dia (12h UT) de 1 de janeiro de 4713 a.C.**
Ela é especial porque marca o início comum de três grandes ciclos de tempo da antiguidade:
1. Ciclo Solar (28 anos)
2. Ciclo Metônico/Lunar (19 anos)
3. Ciclo da Indicção Romana (15 anos)
A última vez que o início destes três ciclos coincidiu foi no ano 4713 a.C.

### B. Por que o dia juliano começa ao meio-dia e não à meia-noite?
Foi definido por astrônomos para que toda a noite de observações ficasse no mesmo dia do calendário juliano. Se começasse à meia-noite, as anotações e medições feitas em uma única noite ficariam divididas em dois dias diferentes.

### C. Quanto tempo passou entre 12/Dez/2008 e 25/mar/2020?
Calculando a Data Juliana para as duas datas:
* $JD(12/12/2008) = 2454813$
* $JD(25/03/2020) = 2458934$
* Dias passados: $2458934 - 2454813 = \mathbf{4121\text{ dias}}$

---

## 3. Explicação das Contas (Algoritmo do Wertz)

O cálculo usa o algoritmo de Fliegel e Van Flandern (1968) baseado em aritmética inteira.
A fórmula divide as contas em 3 partes ($DJ1, DJ2, DJ3$):

$$JD = DJ1 + DJ2 + DJ3$$

Como o MATLAB/Octave calcula divisões normais em ponto flutuante, usamos a função `fix()` para truncar todas as divisões para o menor valor inteiro em direção ao zero, que é a regra da aritmética inteira do FORTRAN antigo.

As equações no código ficaram:
* `DJ1 = K - 32075 + fix((1461 * (I + 4800 + fix((J - 14)/12))) / 4)`
* `DJ2 = fix((367 * (J - 2 - 12 * fix((J - 14)/12))) / 12)`
* `DJ3 = fix((-3 * fix((I + 4900 + fix((J - 14)/12)) / 100)) / 4)`

### Validação dos exemplos:
1. **Para 25/Dez/1981** ($I=1981, J=12, K=25$):
   * $DJ1 = 2444710$
   * $DJ2 = 305$
   * $DJ3 = -51$
   * $JD = 2444964$ (Bate com a resposta esperada)
2. **Para 24/Mar/2020** ($I=2020, J=3, K=24$):
   * $DJ1 = 2458954$
   * $DJ2 = 30$
   * $DJ3 = -51$
   * $JD = 2458933$ (Bate com a resposta esperada)

---

## 4. Comparação com a função `juliandate` do MATLAB

* Se rodar `juliandate(1981, 12, 25)` no MATLAB, o resultado dá **2444963.5**.
* O nosso algoritmo dá **2444964**.

**Por que a diferença?**
A função interna do MATLAB assume que se você não passar a hora, a data refere-se à meia-noite (00h00) do dia civil. Como o dia juliano começa às 12h00 UT, às 00h00 ainda estávamos na metade do dia juliano anterior (por isso dá o .5 a menos).
Para dar igual ao nosso algoritmo, temos que chamar o MATLAB passando as 12 horas: `juliandate(1981, 12, 25, 12, 0, 0)`.

---

## 5. Data Juliana Modificada (DJM) e Fração do Dia

A DJM serve para diminuir o tamanho dos números e sincronizar o início do dia com a meia-noite civil ($00h00 UT$).
O cálculo da DJM no instante de 0h do dia ($DJM0$) é:
$$DJM0 = JD - 2400001$$

Para colocar a fração de horas, minutos e segundos, convertemos o horário em dias:
$$\text{frac\_dia} = \frac{h + \frac{min}{60} + \frac{sec}{3600}}{24}$$
E somamos:
$$DJM_{\text{frac\_dia}} = DJM0 + \text{frac\_dia}$$

* **Exemplo de 9h 10m 35s**:
  $\text{frac\_dia} = [9 + 10*(1/60) + 35*(1/3600)] * (1/24) = 0.382349537037037$ dia.
  $DJM_{\text{frac\_dia}} = 44963 + 0.382349537037037 = \mathbf{44963.382349537037}$ (Exatamente igual ao PDF).

---

## 6. Como rodar os arquivos

Para rodar os testes criados:
1. Abra o MATLAB ou Octave.
2. Navegue até a pasta `atividade_aula_2026_07_01`.
3. Digite e execute:
   ```matlab
   run('testar_DiaJuliano.m')
   ```
Isso mostrará a validação de todos os exemplos e o cálculo da Data Juliana e DJM com fração do dia para o instante de **agora (tempo real)**.
