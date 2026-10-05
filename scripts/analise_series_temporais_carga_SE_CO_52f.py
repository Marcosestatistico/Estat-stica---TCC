from pathlib import Path
import sys, json, warnings, math, platform, itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import friedmanchisquare
import statsmodels
from statsmodels.tsa.stattools import adfuller, kpss, acf, pacf
from statsmodels.tsa.seasonal import STL
from statsmodels.tsa.statespace.sarimax import SARIMAX
from statsmodels.tsa.holtwinters import ExponentialSmoothing
from statsmodels.stats.diagnostic import acorr_ljungbox, het_arch
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')
np.random.seed(20260905)
ROOT = Path(r"C:\Temp\tmp_analysis"); INPUT=Path("CARGA_ENERGIA_CONSOLIDADO.xlsx")
OUT=ROOT/'outputs'; FIG=OUT/'figures'; TAB=OUT/'tables'; FIG.mkdir(parents=True, exist_ok=True); TAB.mkdir(parents=True, exist_ok=True)

# ---------- utilitarios ----------
def aicc(llf, k, n):
    aic=-2*llf+2*k
    return aic + (2*k*(k+1))/(n-k-1) if n>k+1 else np.inf

def metrics(y,p):
    y=np.asarray(y,float); p=np.asarray(p,float); e=y-p
    rmse=np.sqrt(np.mean(e**2)); mae=np.mean(np.abs(e))
    mape=np.mean(np.abs(e/y))*100 if np.all(np.abs(y)>1e-8) else np.nan
    smape=np.mean(200*np.abs(e)/(np.abs(y)+np.abs(p)))
    return {'RMSE':rmse,'MAE':mae,'MAPE':mape,'sMAPE':smape,'ME':np.mean(e)}

def folds_expanding(n, initial, horizon=7, step=28):
    out=[]; end=initial
    while end+horizon<=n:
        out.append((np.arange(0,end),np.arange(end,end+horizon)))
        end += step
    return out

def fourier_matrix(dates, period, K, origin):
    t=(pd.DatetimeIndex(dates)-pd.Timestamp(origin)).days.to_numpy(dtype=float)
    cols={}
    for k in range(1,K+1):
        cols[f'sin_{period}_{k}']=np.sin(2*np.pi*k*t/period)
        cols[f'cos_{period}_{k}']=np.cos(2*np.pi*k*t/period)
    return pd.DataFrame(cols,index=pd.DatetimeIndex(dates))

def adf_kpss(y):
    ad=adfuller(y,autolag='AIC'); kp=kpss(y,regression='ct',nlags='auto')
    return {'ADF_stat':ad[0],'ADF_p':ad[1],'KPSS_stat':kp[0],'KPSS_p':kp[1]}

def select_sarima(y):
    # Grade parcimoniosa predefinida, selecionada por AICc apenas no treino.
    candidates=[((1,0,1),(1,0,1,7)),((1,1,1),(1,0,1,7)),((1,0,1),(0,1,1,7))]
    best=None; failures=[]
    for order,sorder in candidates:
        try:
          d=order[1]; D=sorder[1]
          m=SARIMAX(y,order=order,seasonal_order=sorder,trend='c' if d==0 and D==0 else 'n',
                    enforce_stationarity=False,enforce_invertibility=False).fit(disp=False,maxiter=80)
          score=aicc(m.llf,len(m.params),len(y))
          if np.isfinite(score) and (best is None or score<best[0]): best=(score,order,sorder,m)
        except Exception as ex: failures.append(str(ex)[:120])
    if best is None: raise RuntimeError('Nenhum SARIMA convergiu')
    return best, failures

def select_ets(y):
    best=None; failures=[]
    positive=np.all(np.asarray(y)>0)
    for trend in [None,'add']:
      for damped in ([False,True] if trend else [False]):
       for seasonal in ['add'] + (['mul'] if positive else []):
        try:
          m=ExponentialSmoothing(y,trend=trend,damped_trend=damped,seasonal=seasonal,seasonal_periods=7,
                                 initialization_method='estimated').fit(optimized=True,use_brute=False)
          k=len(m.params); n=len(y); sse=max(m.sse,1e-12)
          ll=-n/2*(np.log(2*np.pi)+1+np.log(sse/n)); score=aicc(ll,k,n)
          if best is None or score<best[0]: best=(score,(trend,damped,seasonal),m)
        except Exception as ex: failures.append(str(ex)[:120])
    if best is None: raise RuntimeError('Nenhum ETS convergiu')
    return best, failures

def select_dhr(y, dates):
    # DHR = regressao harmonica dinamica (Fourier semanal+anual) com erros ARMA/SARIMAX.
    # K anual e ordens de erro selecionados por AICc apenas no treino.
    origin=dates[0]; best=None; failures=[]
    Xw=fourier_matrix(dates,7,3,origin)
    for Ka in [2,4]:
      X=pd.concat([Xw,fourier_matrix(dates,365.25,Ka,origin)],axis=1)
      for p,q in [(0,0),(1,1)]:
       try:
        m=SARIMAX(y,exog=X,order=(p,0,q),seasonal_order=(0,0,0,0),trend='c',
                  enforce_stationarity=False,enforce_invertibility=False).fit(disp=False,maxiter=150)
        score=aicc(m.llf,len(m.params),len(y))
        if np.isfinite(score) and (best is None or score<best[0]): best=(score,Ka,(p,0,q),origin,m)
       except Exception as ex: failures.append(str(ex)[:120])
    if best is None: raise RuntimeError('Nenhum DHR convergiu')
    return best,failures

def dm_test(loss1,loss2,h=7):
    d=np.asarray(loss1)-np.asarray(loss2); n=len(d); mean=d.mean(); z=d-mean
    lag=min(h-1,n-1)
    gamma0=np.dot(z,z)/n; lrv=gamma0
    for k in range(1,lag+1):
        gamma=np.dot(z[k:],z[:-k])/n
        lrv += 2*(1-k/(lag+1))*gamma
    if lrv<=0: return np.nan,np.nan
    stat=mean/np.sqrt(lrv/n)
    # Harvey-Leybourne-Newbold small-sample correction
    corr=np.sqrt((n+1-2*h+h*(h-1)/n)/n) if n>h else 1
    stat*=corr; p=2*stats.t.sf(abs(stat),df=n-1)
    return stat,p

# ---------- leitura e qualidade ----------
xl=pd.ExcelFile(INPUT,engine='openpyxl'); sheet=xl.sheet_names[0]
raw=pd.read_excel(INPUT,sheet_name=sheet,engine='openpyxl')
raw.columns=[str(c).strip() for c in raw.columns]
cols=raw.columns.tolist()
DATE='din_instante'; Y='val_cargaenergiamwmed'; SUB='nom_subsistema'; ID='id_subsistema'
if not all(c in raw.columns for c in [DATE,Y,SUB]): raise ValueError(f'Colunas esperadas nao encontradas. Encontradas: {cols}')
sub_values=raw[SUB].dropna().astype(str).unique().tolist(); filter_value='Sudeste/Centro-Oeste'
if filter_value not in sub_values: raise ValueError(f'Valor exato nao encontrado. Valores: {sub_values}')
df=raw.loc[raw[SUB].astype(str).eq(filter_value)].copy()
# As datas no arquivo sao mm/dd/yyyy, inferido de valores inequivocos como 01/13/2025.
df[DATE]=pd.to_datetime(df[DATE],format='%m/%d/%Y',errors='raise')
df[Y]=pd.to_numeric(df[Y],errors='coerce'); df=df.sort_values(DATE)
duplicates=df.duplicated(DATE,keep=False); full=pd.date_range(df[DATE].min(),df[DATE].max(),freq='D'); missing_dates=full.difference(df[DATE])
quality={'arquivo':INPUT.name,'formato':'xlsx','planilha':sheet,'dimensoes_raw':list(raw.shape),'colunas':cols,
         'tipos_raw':{c:str(t) for c,t in raw.dtypes.items()},'subsistemas':sub_values,'filtro':filter_value,
         'periodo_inicio':str(df[DATE].min().date()),'periodo_fim':str(df[DATE].max().date()),'n':len(df),
         'duplicatas_data':int(duplicates.sum()),'datas_ausentes':len(missing_dates),'nulos':int(df[Y].isna().sum()),
         'negativos':int((df[Y]<0).sum()),'zeros':int((df[Y]==0).sum())}
if quality['duplicatas_data'] or quality['datas_ausentes'] or quality['nulos']: raise ValueError(f'Qualidade impede modelagem sem decisao: {quality}')
ts=df.set_index(DATE)[Y].asfreq('D')
q1,q3=ts.quantile([.25,.75]); iqr=q3-q1; flags=((ts<q1-1.5*iqr)|(ts>q3+1.5*iqr))
quality['outliers_IQR_sinalizados']=int(flags.sum()); quality['limite_IQR_inf']=q1-1.5*iqr; quality['limite_IQR_sup']=q3+1.5*iqr
pd.DataFrame([quality]).to_json(TAB/'qualidade.json',orient='records',indent=2,force_ascii=False)
raw.head().to_csv(TAB/'cinco_primeiras.csv',index=False); raw.tail().to_csv(TAB/'cinco_ultimas.csv',index=False)

desc=ts.describe(percentiles=[.01,.05,.25,.5,.75,.95,.99]).rename('MWmed').to_frame(); desc.to_csv(TAB/'descritivas.csv')
bydow=ts.groupby(ts.index.dayofweek).agg(['count','mean','std','median']); bydow.index=['segunda','terca','quarta','quinta','sexta','sabado','domingo']; bydow.to_csv(TAB/'por_dia_semana.csv')
bymonth=ts.groupby(ts.index.month).agg(['count','mean','std','median']); bymonth.to_csv(TAB/'por_mes.csv')
station=adf_kpss(ts); pd.DataFrame([station]).to_csv(TAB/'estacionariedade.csv',index=False)
stl=STL(ts,period=7,robust=True).fit()

# graficos exploratorios
plt.figure(figsize=(13,5)); plt.plot(ts.index,ts,lw=.7); plt.title('Figura 1. Carga diaria de energia, Sudeste/Centro-Oeste'); plt.ylabel('Carga (MWmed)'); plt.xlabel('Data'); plt.tight_layout(); plt.savefig(FIG/'fig01_serie.png',dpi=180); plt.close()
fig,ax=plt.subplots(4,1,figsize=(13,9),sharex=True); ax[0].plot(ts); ax[0].set_ylabel('MWmed'); ax[0].set_title('Figura 2. Decomposicao STL semanal'); ax[1].plot(stl.trend); ax[1].set_ylabel('Tendencia'); ax[2].plot(stl.seasonal); ax[2].set_ylabel('Sazonal'); ax[3].plot(stl.resid); ax[3].set_ylabel('Residuo'); plt.tight_layout(); plt.savefig(FIG/'fig02_stl.png',dpi=180); plt.close()
fig,ax=plt.subplots(1,2,figsize=(12,4)); lags=60; av=acf(ts,nlags=lags,fft=True); pv=pacf(ts,nlags=lags,method='ywm'); ax[0].stem(range(len(av)),av); ax[0].set_title('Figura 3A. ACF'); ax[0].set_xlabel('Defasagem (dias)'); ax[1].stem(range(len(pv)),pv); ax[1].set_title('Figura 3B. PACF'); ax[1].set_xlabel('Defasagem (dias)'); plt.tight_layout(); plt.savefig(FIG/'fig03_acf_pacf.png',dpi=180); plt.close()
fig,ax=plt.subplots(1,2,figsize=(12,4)); ax[0].hist(ts,bins=35); ax[0].set_title('Figura 4A. Distribuicao'); ax[0].set_xlabel('Carga (MWmed)'); ax[1].boxplot([ts[ts.index.dayofweek==i] for i in range(7)],tick_labels=['Seg','Ter','Qua','Qui','Sex','Sab','Dom']); ax[1].set_title('Figura 4B. Boxplot por dia'); ax[1].set_ylabel('MWmed'); plt.tight_layout(); plt.savefig(FIG/'fig04_distribuicao_boxplot.png',dpi=180); plt.close()

# ---------- validacao ----------
h=7; step=7; nfolds=4; initial=len(ts) - 365; folds=folds_expanding(len(ts),initial,h,step)
fold_rows=[]; pred_rows=[]; spec_rows=[]; diag_rows=[]; failures=[]
for fi,(tri,tei) in enumerate(folds,1):
    train=ts.iloc[tri]; test=ts.iloc[tei]
    assert train.index.max()<test.index.min() and len(test)==h
    fold_rows.append({'fold':fi,'treino_inicio':train.index.min(),'treino_fim':train.index.max(),'teste_inicio':test.index.min(),'teste_fim':test.index.max(),'n_treino':len(train),'n_teste':len(test)})
    # benchmark
    bp=np.array([ts.loc[d-pd.Timedelta(days=7)] for d in test.index])
    for d,a,p in zip(test.index,test,bp): pred_rows.append({'fold':fi,'data':d,'modelo':'SNaive7','real':a,'prev':p,'lo80':np.nan,'hi80':np.nan,'lo95':np.nan,'hi95':np.nan})
    # SARIMA
    try:
      b,fail=select_sarima(train); score,order,sorder,m=b; pr=m.get_forecast(h); mean=np.asarray(pr.predicted_mean); ci80=np.asarray(pr.conf_int(alpha=.20)); ci95=np.asarray(pr.conf_int(alpha=.05))
      spec_rows.append({'fold':fi,'modelo':'SARIMA','especificacao':f'SARIMA{order}x{sorder}','AICc':score})
      r=np.asarray(m.resid)[max(order[1]+7*sorder[1],1):]; lb=acorr_ljungbox(r,lags=[14],return_df=True).iloc[0]
      diag_rows.append({'fold':fi,'modelo':'SARIMA','media_residuo':r.mean(),'dp_residuo':r.std(ddof=1),'LjungBox_Q14':lb.lb_stat,'LjungBox_p14':lb.lb_pvalue})
      for j,(d,a,p) in enumerate(zip(test.index,test,mean)): pred_rows.append({'fold':fi,'data':d,'modelo':'SARIMA','real':a,'prev':p,'lo80':ci80[j,0],'hi80':ci80[j,1],'lo95':ci95[j,0],'hi95':ci95[j,1]})
      failures += [{'fold':fi,'modelo':'SARIMA','erro':x} for x in fail[:2]]
    except Exception as ex: failures.append({'fold':fi,'modelo':'SARIMA','erro':str(ex)})
    # ETS
    try:
      b,fail=select_ets(train); score,conf,m=b; mean=np.asarray(m.forecast(h)); r=np.asarray(m.resid); sig=np.std(r,ddof=1)
      spec_rows.append({'fold':fi,'modelo':'ETS','especificacao':f'ETS(A,{conf[0] or "N"},{conf[2].upper()}; damped={conf[1]})','AICc':score})
      lb=acorr_ljungbox(r,lags=[14],return_df=True).iloc[0]; diag_rows.append({'fold':fi,'modelo':'ETS','media_residuo':r.mean(),'dp_residuo':sig,'LjungBox_Q14':lb.lb_stat,'LjungBox_p14':lb.lb_pvalue})
      for j,(d,a,p) in enumerate(zip(test.index,test,mean)):
        pred_rows.append({'fold':fi,'data':d,'modelo':'ETS','real':a,'prev':p,'lo80':p-stats.norm.ppf(.9)*sig,'hi80':p+stats.norm.ppf(.9)*sig,'lo95':p-1.96*sig,'hi95':p+1.96*sig})
      failures += [{'fold':fi,'modelo':'ETS','erro':x} for x in fail[:2]]
    except Exception as ex: failures.append({'fold':fi,'modelo':'ETS','erro':str(ex)})
    # DHR Fourier + ARMA errors
    try:
      b,fail=select_dhr(train,train.index); score,Ka,order,origin,m=b
      Xtest=pd.concat([fourier_matrix(test.index,7,3,origin),fourier_matrix(test.index,365.25,Ka,origin)],axis=1)
      pr=m.get_forecast(h,exog=Xtest); mean=np.asarray(pr.predicted_mean); ci80=np.asarray(pr.conf_int(alpha=.20)); ci95=np.asarray(pr.conf_int(alpha=.05))
      spec_rows.append({'fold':fi,'modelo':'DHR-Fourier','especificacao':f'Fourier semanal K=3 + anual K={Ka} + ARMA{(order[0],order[2])}','AICc':score})
      r=np.asarray(m.resid); lb=acorr_ljungbox(r,lags=[14],return_df=True).iloc[0]; diag_rows.append({'fold':fi,'modelo':'DHR-Fourier','media_residuo':r.mean(),'dp_residuo':r.std(ddof=1),'LjungBox_Q14':lb.lb_stat,'LjungBox_p14':lb.lb_pvalue})
      for j,(d,a,p) in enumerate(zip(test.index,test,mean)): pred_rows.append({'fold':fi,'data':d,'modelo':'DHR-Fourier','real':a,'prev':p,'lo80':ci80[j,0],'hi80':ci80[j,1],'lo95':ci95[j,0],'hi95':ci95[j,1]})
      failures += [{'fold':fi,'modelo':'DHR-Fourier','erro':x} for x in fail[:2]]
    except Exception as ex: failures.append({'fold':fi,'modelo':'DHR-Fourier','erro':str(ex)})
    print('fold',fi,'ok')

fold_df=pd.DataFrame(fold_rows); pred=pd.DataFrame(pred_rows); specs=pd.DataFrame(spec_rows); diags=pd.DataFrame(diag_rows)
fold_df.to_csv(TAB/'folds.csv',index=False); pred.to_csv(TAB/'previsoes.csv',index=False); specs.to_csv(TAB/'especificacoes.csv',index=False); diags.to_csv(TAB/'diagnosticos.csv',index=False); pd.DataFrame(failures).to_csv(TAB/'falhas.csv',index=False)
expected = len(folds) * h; counts = pred.groupby('modelo').size(); assert all(counts == expected), counts

fold_metrics=[]
for (f,m),g in pred.groupby(['fold','modelo']): fold_metrics.append({'fold':f,'modelo':m,**metrics(g.real,g.prev)})
fm=pd.DataFrame(fold_metrics); fm.to_csv(TAB/'metricas_por_fold.csv',index=False)
summary=[]
for m,g in fm.groupby('modelo'):
  gp=pred[pred.modelo==m]; ag=metrics(gp.real,gp.prev)
  for met in ['RMSE','MAE','MAPE','sMAPE']:
    vals=g[met]; ci=stats.t.interval(.95,len(vals)-1,loc=vals.mean(),scale=stats.sem(vals))
    summary.append({'modelo':m,'metrica':met,'media_folds':vals.mean(),'mediana_folds':vals.median(),'dp_folds':vals.std(ddof=1),'IC95_inf':ci[0],'IC95_sup':ci[1],'agregada_todas_previsoes':ag[met]})
sumdf=pd.DataFrame(summary); sumdf.to_csv(TAB/'resumo_metricas.csv',index=False)
rank=sumdf[sumdf.metrica=='RMSE'].sort_values('media_folds').copy(); rank['ranking_RMSE']=range(1,len(rank)+1); rank.to_csv(TAB/'ranking.csv',index=False)
# Cobertura e largura
ints=[]
for m,g in pred[pred.modelo!='SNaive7'].groupby('modelo'):
  for lev in [80,95]:
    lo=g[f'lo{lev}']; hi=g[f'hi{lev}']; ints.append({'modelo':m,'nivel':lev,'cobertura':np.mean((g.real>=lo)&(g.real<=hi))*100,'largura_media':np.mean(hi-lo)})
intdf=pd.DataFrame(ints); intdf.to_csv(TAB/'intervalos.csv',index=False)
# DM squared error, mesma data, Holm
wide=pred.pivot(index='data',columns='modelo',values=['real','prev']); models=['SARIMA','ETS','DHR-Fourier']; dm=[]
for a,b in itertools.combinations(models,2):
 y=wide['real'][a].to_numpy(); l1=(y-wide['prev'][a].to_numpy())**2; l2=(y-wide['prev'][b].to_numpy())**2
 st,p=dm_test(l1,l2,h); dm.append({'modelo_1':a,'modelo_2':b,'DM_stat':st,'p_bruto':p,'perda_media_1_menos_2':np.mean(l1-l2)})
padj=multipletests([x['p_bruto'] for x in dm],alpha=.05,method='holm')[1]
for x,p in zip(dm,padj): x['p_Holm']=p; x['significativo_5pct']=bool(p<.05)
dmdf=pd.DataFrame(dm); dmdf.to_csv(TAB/'diebold_mariano.csv',index=False)
# Friedman em RMSE dos folds (complementar)
piv=fm[fm.modelo.isin(models)].pivot(index='fold',columns='modelo',values='RMSE'); fr=friedmanchisquare(*[piv[m] for m in models]); pd.DataFrame([{'estatistica':fr.statistic,'p':fr.pvalue,'n_folds':len(piv)}]).to_csv(TAB/'friedman.csv',index=False)

# graficos de resultados
plt.figure(figsize=(13,5));
for m,g in pred.groupby('modelo'): plt.plot(g.data,g.prev,label=m,lw=1)
plt.scatter(pred[pred.modelo=='SARIMA'].data,pred[pred.modelo=='SARIMA'].real,s=10,c='black',label='Real'); plt.ylabel('Carga (MWmed)'); plt.xlabel('Data'); plt.title('Figura 5. Previsoes fora da amostra'); plt.legend(ncol=5,fontsize=8); plt.tight_layout(); plt.savefig(FIG/'fig05_real_previsto.png',dpi=180); plt.close()
err=[(pred[pred.modelo==m].real-pred[pred.modelo==m].prev).to_numpy() for m in ['SARIMA','ETS','DHR-Fourier','SNaive7']]; plt.figure(figsize=(9,5)); plt.boxplot(err,tick_labels=['SARIMA','ETS','DHR-Fourier','SNaive7'],showfliers=True); plt.ylabel('Erro real - previsto (MWmed)'); plt.title('Figura 6. Distribuicao dos erros fora da amostra'); plt.tight_layout(); plt.savefig(FIG/'fig06_erros.png',dpi=180); plt.close()
plt.figure(figsize=(10,5));
for m,g in fm.groupby('modelo'): plt.plot(g.fold,g.RMSE,marker='o',label=m)
plt.ylabel('RMSE (MWmed)'); plt.xlabel('Fold'); plt.title('Figura 7. RMSE por fold'); plt.legend(); plt.tight_layout(); plt.savefig(FIG/'fig07_rmse_folds.png',dpi=180); plt.close()
# final fit and 7-day forecast for illustrative final model intervals
winner=rank.iloc[0].modelo; train=ts
if winner=='SARIMA':
 b,_=select_sarima(train); _,order,sorder,m=b; fut=pd.date_range(ts.index[-1]+pd.Timedelta(days=1),periods=7); pr=m.get_forecast(7); mean=np.asarray(pr.predicted_mean); c80=np.asarray(pr.conf_int(alpha=.2)); c95=np.asarray(pr.conf_int(alpha=.05)); spec=f'SARIMA{order}x{sorder}'
elif winner=='ETS':
 b,_=select_ets(train); _,conf,m=b; fut=pd.date_range(ts.index[-1]+pd.Timedelta(days=1),periods=7); mean=np.asarray(m.forecast(7)); sig=np.std(m.resid,ddof=1); c80=np.c_[mean-stats.norm.ppf(.9)*sig,mean+stats.norm.ppf(.9)*sig]; c95=np.c_[mean-1.96*sig,mean+1.96*sig]; spec=f'ETS(A,{conf[0] or "N"},{conf[2]})'
else:
 b,_=select_dhr(train,train.index); _,Ka,order,origin,m=b; fut=pd.date_range(ts.index[-1]+pd.Timedelta(days=1),periods=7); X=pd.concat([fourier_matrix(fut,7,3,origin),fourier_matrix(fut,365.25,Ka,origin)],axis=1); pr=m.get_forecast(7,exog=X); mean=np.asarray(pr.predicted_mean); c80=np.asarray(pr.conf_int(alpha=.2)); c95=np.asarray(pr.conf_int(alpha=.05)); spec=f'DHR K7=3 K365={Ka} ARMA{(order[0],order[2])}'
final=pd.DataFrame({'data':fut,'prev':mean,'lo80':c80[:,0],'hi80':c80[:,1],'lo95':c95[:,0],'hi95':c95[:,1]}); final.to_csv(TAB/'previsao_final.csv',index=False)
plt.figure(figsize=(12,5)); recent=ts.tail(90); plt.plot(recent.index,recent,label='Observado'); plt.plot(fut,mean,label=f'Previsao {winner}'); plt.fill_between(fut,c95[:,0],c95[:,1],alpha=.15,label='IP 95%'); plt.fill_between(fut,c80[:,0],c80[:,1],alpha=.25,label='IP 80%'); plt.ylabel('Carga (MWmed)'); plt.xlabel('Data'); plt.title('Figura 8. Previsao final de 7 dias'); plt.legend(); plt.tight_layout(); plt.savefig(FIG/'fig08_previsao_final.png',dpi=180); plt.close()

# maiores erros
pred['erro_abs']=(pred.real-pred.prev).abs(); pred.sort_values('erro_abs',ascending=False).head(20).to_csv(TAB/'maiores_erros.csv',index=False)
results={'quality':quality,'stationarity':station,'winner':winner,'winner_spec_final':spec,
         'ranking':rank.to_dict('records'),'dm':dm,'friedman':{'stat':fr.statistic,'p':fr.pvalue},
         'intervals':ints,'counts':counts.to_dict(),'python':sys.version,'pandas':pd.__version__,'numpy':np.__version__,'scipy':stats.__version__ if hasattr(stats,'__version__') else __import__('scipy').__version__,'statsmodels':statsmodels.__version__}
with open(OUT/'results.json','w',encoding='utf-8') as f: json.dump(results,f,ensure_ascii=False,indent=2,default=str)
print(json.dumps(results,ensure_ascii=False,indent=2,default=str))
