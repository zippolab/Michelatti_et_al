import graph_tool.all as gt
import numpy as np
import pandas as pd


#		1. COMPUTING BASIC STATISTICS OVER NASCENT CLUSTERS NETWORKS


ntw_cat=['vertex_n','edges_n','chr','peak','cl1','cl2','cl3','cl4','cl5','cl6','not_anno']

for conn in ['cc','cov']:
    for thr in [0,1,2,3,4,5,6,7,8,9]:
        for cl in tot_cl['cl'].cat.categories:
            DF=pd.DataFrame(columns=ntw_cat)
            df_cl=tot_cl[tot_cl['cl']==cl]
            for gn in df_cl[8].tolist():
                try:
                    g=gt.load_graph(f'graphs/{cl}/{cl}_{gn}_g_cluster0_{thr}perc_tot_{conn}_linked.xml.gz')
                    ref=pd.read_table(f'dataframes/{cl}/{cl}_{gn}_reg_cluster0_{thr}perc_tot_{conn}_linked.bed')
                    cluster=[]
                    for cat in 'cl1','cl2','cl3','cl4','cl5','cl6','0':
                        cluster.append(ref[ref['nascent_cl']==cat].shape[0])
                    stats=[g.num_vertices(),g.num_edges(),ref['chr'].tolist()[0].split('chr')[1],sum(ref['peak'])/ref.shape[0]]
                    stats=stats+cluster
                    df=pd.DataFrame(stats,columns=[gn],index=ntw_cat).T
                    df['thr']=[thr for x in range(df.shape[0])]
        	    df['gene']=df.index
        	    df['cluster']=[cl for x in range(df.shape[0])]
                    DF=pd.concat([DF,df])
                except FileNotFoundError:
                    df=pd.DataFrame([0 for x in range(len(ntw_cat))],columns=[gn],index=ntw_cat).T
                    df['thr']=[thr for x in range(df.shape[0])]
        	    df['gene']=df.index
        	    df['cluster']=[cl for x in range(df.shape[0])]
                    DF=pd.concat([DF,df])
            DF.to_csv(f'analysis/{conn}_filtered/{cl}/{cl}_{thr}perc_{conn}_basic_measures_linked.bed',sep='\t')


#		2. EVALUATING NASCENT PROMOTERS EROSION FOR EACH NASCENT GROUP


DF_PLOT=pd.DataFrame()
for cl in tot_cl['cl'].cat.categories:
    df_plot=pd.DataFrame(columns=['relative abundance','percentile','nascent cluster'])
    for thr in [0,1,2,3,4,5,6,7,8,9]:
        df_cl=DF[(DF['cluster']==cl)&(DF['thr']==thr)]
        df_cl=df_cl.iloc[:,7:-1]
        df_cl=df_cl/sum(sum(df_cl))
        df_cl_plot=pd.DataFrame([[x for x in sum(df_cl)],[thr for x in range(6)],['cl1', 'cl2', 'cl3', 'cl6', 'mem', 'res']],index=['relative abundance','percentile','nascent cluster']).T
        df_plot=pd.concat([df_plot,df_cl_plot])
    df_plot.index=[x for x in range(df_plot.shape[0])]
    df_plot=df_plot[df_plot['nascent cluster']==cl]
    DF_PLOT=pd.concat([DF_PLOT,df_plot])
DF_PLOT.index=[x for x in range(DF_PLOT.shape[0])]

#		3. FOCUS ON PROMOTER EROSION IN MEMORY AND RESPONSIVE ONLY CLUSTERS


df_plot_mem=DF_PLOT[DF_PLOT['nascent cluster']=='mem']
df_plot_res=DF_PLOT[DF_PLOT['nascent cluster']=='res']

df_plot_mem['relative abundance']=df_plot_mem['relative abundance']/max(df_plot_mem['relative abundance'])
df_plot_res['relative abundance']=df_plot_res['relative abundance']/max(df_plot_res['relative abundance'])

DF_PLOT=pd.concat([df_plot_mem,df_plot_res])
