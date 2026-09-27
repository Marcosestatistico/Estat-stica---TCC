# Comparação de modelos de séries temporais para previsão da carga diária de energia elétrica

## Tema geral

Previsão estatística de demanda energética.

## Problema de pesquisa

Qual modelo produz melhor desempenho preditivo para a carga diária de energia do subsistema Sudeste/Centro-Oeste: **ARIMA**, **ETS** ou **DHR**?

## Objetivo geral

Comparar o desempenho de modelos de séries temporais para previsão da carga diária de energia.

## Objetivos específicos

1. Verificar se modelos que incorporam tendência e sazonalidade apresentam melhor desempenho de predição;
2. Verificar qual modelo produz menor erro preditivo para a carga diária de energia, utilizando como métricas **RMSE**, **MAE** e **MAPE**;
3. Verificar se há diferenças estatisticamente significativas no desempenho preditivo dos modelos;
4. Fazer recomendações para previsão da carga diária de energia.

## Hipóteses

### Hipótese científica

Modelos que incorporam adequadamente tendência e sazonalidade apresentam melhor desempenho preditivo para a carga diária de energia.

### Hipóteses estatísticas

- **H0:** Não existem diferenças estatisticamente significativas no desempenho preditivo dos modelos.
- **H1:** Pelo menos um modelo apresenta desempenho preditivo estatisticamente diferente dos demais.

## Delimitação

- Subsistema Sudeste/Centro-Oeste do Sistema Interligado Nacional (SIN);
- Dados diários;
- Últimos 5 anos.

## Fonte dos dados

- **Instituição:** Operador Nacional do Sistema Elétrico (ONS);
- **Portal:** ONS Dados Abertos;
- **Conjunto de dados:** Carga de Energia Diária;
- **Acesso:** [Carga de Energia Diária - ONS Dados Abertos](https://dados.ons.org.br/dataset/carga-energia).

> O conjunto disponibiliza dados de carga por subsistema em base diária, medidos em MWmed. Os dados estão sujeitos a um processo recorrente de consistência e podem ser atualizados após a publicação.

## Modelos avaliados

- ARIMA;
- ETS;
- DHR.

## Métricas de avaliação

- RMSE;
- MAE;
- MAPE.

## Palavras-chave

`séries temporais` · `previsão de demanda` · `energia elétrica` · `ARIMA` · `forecasting` · `electric load forecasting` · `ETS` · `forecast evaluation`

## Licença e atribuição dos dados

Os dados utilizados são provenientes do ONS Dados Abertos. Ao reutilizá-los, devem ser observados os termos e as condições informados no portal, incluindo o crédito apropriado ao ONS e a indicação de eventuais alterações realizadas nos dados.
