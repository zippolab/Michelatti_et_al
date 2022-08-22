import pandas as pd
import numpy as np


conn_MD=pd.read_table('graph_anno/MD_tot_prom_anno.bed')

### conn_MD is a dataframe with 5 columns for each HiChIP bins associated with a promoter, describing 3 features associated with the promoter: degree, coverage and connectivity

###bin 	degree 	coverage 	connectivity	symbol





#		1. Generation of the dataframe for kde plot

def kde_percentiles(cell,df,counts,feature,thr,scale):
    cell=cell
    measure=feature
    ranges=[[np.percentile(np.array(df[f'{measure}']),0),np.percentile(np.array(df[f'{measure}']),10)],[np.percentile(np.array(df[f'{measure}']),10),np.percentile(np.array(df[f'{measure}']),20)],[np.percentile(np.array(df[f'{measure}']),20),np.percentile(np.array(df[f'{measure}']),30)],[np.percentile(np.array(df[f'{measure}']),30),np.percentile(np.array(df[f'{measure}']),40)],[np.percentile(np.array(df[f'{measure}']),40),np.percentile(np.array(df[f'{measure}']),50)],[np.percentile(np.array(df[f'{measure}']),50),np.percentile(np.array(df[f'{measure}']),60)],[np.percentile(np.array(df[f'{measure}']),60),np.percentile(np.array(df[f'{measure}']),70)],[np.percentile(np.array(df[f'{measure}']),70),np.percentile(np.array(df[f'{measure}']),80)],[np.percentile(np.array(df[f'{measure}']),80),np.percentile(np.array(df[f'{measure}']),90)],[np.percentile(np.array(df[f'{measure}']),90),np.percentile(np.array(df[f'{measure}']),100)]]

    expression_tot=[]
    for s,e in ranges:
        df_i=df[(df[f'{measure}']>=s)&(df[f'{measure}']<=e)]
        expression_tot.append(counts[counts.index.isin(df_i['symbol'].tolist())][f'{cell}'].tolist())        

    df_sns=pd.DataFrame(columns=['counts','percentile'])
    for i in range(len(expression_tot)):
        exp=expression_tot[i]
        perc5=np.percentile(np.array(exp), thr)
        perc95=np.percentile(np.array(exp),100-thr)
        exp_f=[x for x in exp if x>perc5]
        exp_f=[x for x in exp_f if x<perc95]
        exp=exp_f
        connections=[int((i+1)*10) for x in range(len(exp))]
        df_sns_i=pd.DataFrame([exp,connections],index=['counts','percentile']).T
        df_sns=pd.concat([df_sns,df_sns_i])
    df_sns['counts']=df_sns['counts'].astype(float)
    df_sns_f=df_sns
    return df_sns_f
    
df_kde = kde_percentiles('MD',conn_MD,rna_counts,'connectivity',5,'log')
