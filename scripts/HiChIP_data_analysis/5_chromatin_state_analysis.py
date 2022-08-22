import graph_tool.all as gt
import numpy as np
import pandas as pd



epi_cat=['super_enhancer','broad_peaks','promoter','enhancer','insulator','polycomb-repressed','heterochromatin','other']  # epigenetic labels of HiChIP bins




#		1. COUNTING CHROMATIN STATES LABELLED BINS FOR EACH NASCENT CLUSTER 


for conn in ['cc','cov']:                   # same process for connectivity and coverage percentiles
    for thr in [0,1,2,3,4,5,6,7,8,9]:
        for cl in tot_cl['cl'].cat.categories:
            DF=pd.DataFrame(columns=epi_cat)
            df_cl=tot_cl[tot_cl['cl']==cl]
            for gn in df_cl[8].tolist():
                try:
                    g=gt.load_graph(f'graphs/{cl}/{cl}_{gn}_g_cluster0_{thr}perc_tot_{conn}_linked.xml.gz')
                    ref=pd.read_table(f'dataframes/{cl}/{cl}_{gn}_reg_cluster0_{thr}perc_tot_{conn}_linked.bed')
                    stats=[]
                    for cat in epi_cat:
                        try:
                            count=ref['annotation'].value_counts()[f'{cat}']
                        except KeyError:
                            count=0
                        stats.append(count)
                    df=pd.DataFrame(stats,columns=[gn],index=epi_cat).T
                    DF=pd.concat([DF,df])
                except FileNotFoundError:
                    print(cl,gn)
                    df=pd.DataFrame([0 for x in range(len(epi_cat))],columns=[gn],index=epi_cat).T
                    DF=pd.concat([DF,df])
                DF.to_csv(f'analysis_{conn}_filtered/{cl}/{cl}_{thr}perc_{conn}_chromatin_states_linked.bed',sep='\t')
                
                
#		2. EVALUATE THE FRACTION OF EACH CHROMATIN STATE LABEL IN CONSIDERING ALL CLUSTERS OF EACH NASCENT GROUP
                
conn='cov' #here for coverage, but also for connectivity

DF=pd.DataFrame(columns=epi_cat)

for cl in tot_cl['cl'].cat.categories:
    for thr in [0,1,2,3,4,5,6,7,8,9]:
        df=pd.read_table(f'analysis/{conn}_filtered/{cl}/{cl}_{thr}perc_{conn}_chromatin_states_tot_linked.bed', index_col=0)
        mean_epi=df.sum()/sum(sum(df))
        df_mean=pd.DataFrame(mean_epi,index=mean_epi.index, columns=[f'{cl}_{thr}']).T
        DF=pd.concat([DF,df_mean])
DF=DF.astype(float)


#		3. RESILIENCE OF CHROMATIN STATE ANNOTATION

DF_0=DF[DF.index.str.endswith('0')] #annotations without filtering according to connectivity or coverage
DF_9=DF[DF.index.str.endswith('9')] #annotations, using a filtering thrshold equal to the 90th percentile fo connectivity or coverage
DF_9.index=['cl1','cl2','cl3','cl4','cl5','cl6']
DF_ratio=DF_9/DF_0    # estimating resilience of chromatin annotation

