# Tempos entre falhas de ventiladores pulmonares — dados e scripts

Dados e código do artigo *Weibull Reliability Modelling of Mechanical Ventilators in a Regional Reference Hospital: Failure Characterisation and Transition from Corrective to Reliability-Centred Maintenance* (submetido a *Socio-Economic Planning Sciences*).

*English version below.*

## Português

### O que são os dados
`data/tbf_ventiladores.csv` contém os **100 tempos entre falhas (TBF)**, em **dias**, de **18 ventiladores pulmonares** de um hospital de referência da Macrorregião de Saúde do Vale do Jequitinhonha (Minas Gerais, Brasil). Os dados vêm dos históricos de manutenção corretiva de **quatro anos consecutivos** de operação.

- Cada valor é o intervalo entre duas falhas consecutivas **do mesmo equipamento**. Os 100 intervalos dos 18 equipamentos foram agrupados numa só amostra.
- Todos os intervalos são completos (terminam numa falha registrada); não há observações censuradas.
- O arquivo não identifica qual equipamento gerou cada intervalo, e a coluna `ordem` não indica ordem cronológica.
- Os registros não trazem dados de pacientes. O hospital, o fabricante e o modelo dos equipamentos não são identificados.

### Dicionário de variáveis

| Coluna | Tipo | Descrição |
|---|---|---|
| `ordem` | inteiro (1–100) | Número sequencial do registro |
| `tbf_dias` | decimal | Tempo entre falhas consecutivas do mesmo ventilador, em dias |

### Scripts

| Arquivo | Conteúdo |
|---|---|
| `analise_confiabilidade_ventiladores.R` | Análise principal: ajuste MLE de Weibull, lognormal, exponencial e normal; χ², Kolmogorov–Smirnov, AIC/BIC; IC bootstrap (B = 5.000); R(t), F(t), h(t), H(t); percentis t10/t50/t90; função de custo de ciclo g(T) do modelo de substituição por idade e intervalo ótimo T\*; comparação de estratégias. Gera as figuras em `figs/` (rótulos em português, valores em R$). |
| `fig4_custo_economia_usd.py` | Refaz a Figura 4 (a e b) do artigo em dólares (USD 1,00 = BRL 5,11), com rótulos em inglês. |

O script R traz o vetor de dados embutido (idêntico ao CSV) e roda sem o CSV.

### Como reproduzir
```bash
Rscript analise_confiabilidade_ventiladores.R   # grava figs/*.png e figs/*.pdf (~1 min)
python3 fig4_custo_economia_usd.py              # grava figs/fig4a-*.png/pdf e figs/fig4b-*.png/pdf
```
A semente é fixa (`set.seed(2024)`). Resultados principais: β̂ = 1,445; η̂ = 21,53 dias; MTTF = 19,53 dias; T\* = 14,7 dias.

### Ambiente usado
R 4.6.1 (Ubuntu): fitdistrplus 1.2-6, ggplot2 4.0.3, dplyr 1.2.1, tidyr 1.3.2, patchwork 1.3.2, scales 1.4.0, boot 1.3-32, nortest 1.0-4.
Python 3.14: numpy 2.5.3, matplotlib 3.11.2.
O cabeçalho do script R lista as versões da execução original (R 4.4.1).

```r
install.packages(c("fitdistrplus","ggplot2","dplyr","tidyr","patchwork","scales","boot","nortest"))
```
```bash
pip install numpy matplotlib
```

### Licenças
Dados (`data/tbf_ventiladores.csv`): **CC BY 4.0**. Código (`*.R`, `*.py`): **MIT**. Ver `LICENSE`.

### Como citar
Abreu, C. M.; Moreira, L. C.; Carneiro, G. H. F.; Andrade, P. C. R. (2026). *Tempos entre falhas de ventiladores pulmonares — dados e scripts* [conjunto de dados]. GitHub. https://github.com/guihenrique13/ventilator-reliability-weibull

---

## English

### What the data are
`data/tbf_ventiladores.csv` holds the **100 times between failures (TBF)**, in **days**, of **18 mechanical ventilators** at a reference hospital of the Jequitinhonha Valley Health Macro-region (Minas Gerais, Brazil). The data come from corrective maintenance histories covering **four consecutive years** of operation.

- Each value is the interval between two consecutive failures **of the same device**. The 100 intervals from the 18 devices were pooled into one sample.
- All intervals are complete (each ends in a recorded failure); there are no censored observations.
- The file does not identify which device produced each interval, and the `ordem` column does not indicate chronological order.
- The records contain no patient data. The hospital and the equipment make and model are not identified.

### Variable dictionary

| Column | Type | Description |
|---|---|---|
| `ordem` | integer (1–100) | Sequential record number |
| `tbf_dias` | decimal | Time between consecutive failures of the same ventilator, in days |

### Scripts

| File | Content |
|---|---|
| `analise_confiabilidade_ventiladores.R` | Main analysis: MLE fit of Weibull, lognormal, exponential and normal; χ², Kolmogorov–Smirnov, AIC/BIC; bootstrap CIs (B = 5,000); R(t), F(t), h(t), H(t); t10/t50/t90 percentiles; age-replacement cycle cost g(T) and optimal interval T\*; strategy comparison. Writes figures to `figs/` (Portuguese labels, values in BRL). |
| `fig4_custo_economia_usd.py` | Redraws Figure 4 (a and b) of the paper in US dollars (USD 1.00 = BRL 5.11), with English labels. |

The R script embeds the data vector (identical to the CSV) and runs without the CSV.

### How to reproduce
```bash
Rscript analise_confiabilidade_ventiladores.R   # writes figs/*.png and figs/*.pdf (~1 min)
python3 fig4_custo_economia_usd.py              # writes figs/fig4a-* and figs/fig4b-*
```
The seed is fixed (`set.seed(2024)`). Main results: β̂ = 1.445; η̂ = 21.53 days; MTTF = 19.53 days; T\* = 14.7 days.

### Environment
R 4.6.1 (Ubuntu): fitdistrplus 1.2-6, ggplot2 4.0.3, dplyr 1.2.1, tidyr 1.3.2, patchwork 1.3.2, scales 1.4.0, boot 1.3-32, nortest 1.0-4.
Python 3.14: numpy 2.5.3, matplotlib 3.11.2.

### Licences
Data (`data/tbf_ventiladores.csv`): **CC BY 4.0**. Code (`*.R`, `*.py`): **MIT**. See `LICENSE`.

### How to cite
Abreu, C. M.; Moreira, L. C.; Carneiro, G. H. F.; Andrade, P. C. R. (2026). *Times between failures of mechanical ventilators — data and scripts* [dataset]. GitHub. https://github.com/guihenrique13/ventilator-reliability-weibull
