import graph_tool.all as gt
import numpy as np
import pandas as pd
%pylab
%matplotlib inline


### tot_cl = dataframe of 11 columns. Example:
###  column_0 (bin_chr): chr10, column_1 (bin_start): 131900000, column_2 (bin_end): 131910000, 
### column_3 (promoter_chr): chr10, column_4 (promoter_start): 131901773, column_5 (promoter_end): 131904373,
### column_6 (ENS_GENE id): ENSG00000233122, column_7 (strand): +, column_8 (GENE symbol): CTAGE7P,
### column_9 (bin): chr10-131900000-131910000, column_10 (nascent_cluster): 'cl1'

### anno = dataframe, comprising each HiChIP bin, associated with a chromatin state (i.e. promoter), gene promoters (if any) and nascent cluster annotation (if any)

### thr = dataframe with percentiles of connectivity and coverage of HiChIP experiment

# 	1. iCD GENERATION OF EACH NASCENT GENE
for cl in tot_cl['cl'].cat.categories:
    print(cl)
    df_cl=tot_cl[tot_cl['cl']==cl]
    for i in range(df_cl.shape[0]):
        df=pd.read_csv(f'home/MD-6-H3K27ac_10kb.interactions_FitHiC_Q0.01.bed', delimiter='\t')
        region=df_cl.iloc[i,9]
        gn=df_cl.iloc[i,8]
        df,reg=prep_df(df,False,[region],True)
        g=get_graph(df,reg)
        #ANNOTATION
        # coverage + peaks
        list_of_cov=[]
        list_of_peaks=[]
        for i in reg['regions'].tolist():
            DF=df[df['start']==i]
            if DF.shape[0]>0:
                list_of_cov.append(DF['coverage1'].tolist()[0])
                list_of_peaks.append(DF['peak1'].tolist()[0])
            else:
                DF=df[df['end']==i]
                list_of_cov.append(DF['coverage2'].tolist()[0])
                list_of_peaks.append(DF['peak2'].tolist()[0])

        reg['coverage']=list_of_cov
        reg['peak']=list_of_peaks
        #chromatin annotation,nascent_gene ,nascent cl
        anno_cl=anno[anno['bin'].isin(reg['regions'].tolist())]
        anno_cl=anno_cl.sort_values('bin')
        reg=reg.sort_values('regions')
        if anno_cl['bin'].tolist()==reg['regions'].tolist():
            reg['annotation']=anno_cl['annotation'].tolist()
            reg['nascent_cl']=anno_cl['nascent_cl'].tolist()
            reg['nascent_gene']=anno_cl['nascent_gene'].tolist()
            reg=reg.sort_index()
            df.to_csv(f"dataframes/{cl}/{cl}_{gn}_df_tot.bed",sep='\t',index=False)
            reg.to_csv(f"dataframes/{cl}/{cl}_{gn}_reg_tot.bed",sep='\t',index=False)
            g.save(f"graphs/{cl}/{cl}_{gn}_g_tot.xml.gz")
        else:
            continue
            
            
# 	2. CLUSTER ANALYSIS OF EACH iCD

for cl in tot_cl['cl'].cat.categories:
    print(cl)
    df_cl=tot_cl[tot_cl['cl']==cl]
    for i in range(df_cl.shape[0]):
        region=df_cl.iloc[i,9]
        gn=df_cl.iloc[i,8]
        try:
            try:
                reg=pd.read_table(f"dataframes/{cl}/{cl}_{gn}_reg_tot.bed")
                df=pd.read_table(f"dataframes/{cl}/{cl}_{gn}_df_tot.bed")
                df_linked=link_graph(df,reg)
                g_linked=get_graph(df_linked,reg)
                state_l = gt.minimize_nested_blockmodel_dl(g_linked,state_args=dict(base_type=gt.LayeredBlockState,state_args=dict(ec=g_linked.ep.layer, layers=True)))
                levels = state_l.get_levels()
                #keep only level 0 of nested model
                state_level0=levels[0]
                bl_0=state_level0.get_blocks()  
                # storin of the proper cluster of level 0 (cluster in which the nascent promoter is present)
                # gene vertex code in order to select the proper cluster
                gn_code=reg[reg['regions']==region]['codes'].tolist()[0] 
                cl_code=bl_0.a[gn_code]
                #selection of vertices of cluster of nascent promoter
                list_v=[]
                for i in range(len(bl_0.a)):
                    if bl_0.a[i]==cl_code:
                        list_v.append(i)
                #creation of a mask (true if the corresponding vertex is in the selected cluster, else false)
                mask=[]
                for v in [int(x) for x in g.vertices()]:
                    if v in list_v:
                        mask.append(True)
                    else:
                        mask.append(False)
                # generation of df and reg for the selected cluster
                reg_cluster=reg[reg['codes'].isin(list_v)]
                df_cluster=df[(df['start_codes'].isin(list_v))&(df['end_codes'].isin(list_v))]
                df_cluster.index=[x for x in range(df_cluster.shape[0])]
                reg_cluster.index=[x for x in range(reg_cluster.shape[0])]
                df_cluster.to_csv(f'dataframes/{cl}/{cl}_{gn}_df_cluster0_linked.bed',sep='\t',index=False)
                reg_cluster.to_csv(f'dataframes/{cl}/{cl}_{gn}_reg_cluster0_linked.bed',sep='\t',index=False)
                #generation of a graph comprising only the selected cluster
                g_cluster=get_graph(df_cluster,reg_cluster)
                g_cluster.save(f"graphs/{cl}/{cl}_{gn}_g_cluster0_linked.xml.gz")
            except IndexError:
                print('Ierror',cl,gn)
        except FileNotFoundError: #files whose chromHMM annotations for at least one bin is missing
            print(cl,gn)
            continue


# 	3. CLUSTER FILTERING ACCORDING TO CONNECTIVITY PERCENTILES

flt='connectivity'

for cl in tot_cl['cl'].cat.categories:
    print(cl)
    df_cl=tot_cl[tot_cl['cl']==cl]
    for i in range(df_cl.shape[0]):
        region=df_cl.iloc[i,9]
        gn=df_cl.iloc[i,8]
        try:
            reg=pd.read_table(f"dataframes/{cl}/{cl}_{gn}_reg_cluster0_linked.bed")
            df=pd.read_table(f"dataframes/{cl}/{cl}_{gn}_df_cluster0_linked.bed")
            # thresholds definition
            thr_tot=thr.T[f'{flt}_tot'].astype(float32).tolist()            
            for t in range(len(thr_tot)):
                tr=thr_tot[t]
                df_t=df[df['contacts']>=tr]
                if df_t.shape[0]>=1:
                    reg_t=reg[reg['regions'].isin(list(set(df_t['start'])|set(df_t['end'])))]
                    df_t.to_csv(f'dataframes/{cl}/{cl}_{gn}_df_cluster0_{t}perc_tot_{flt}_linked.bed',sep='\t',index=False)
                    reg_t.to_csv(f'dataframes/{cl}/{cl}_{gn}_reg_cluster0_{t}perc_tot_{flt}_linked.bed',sep='\t',index=False)
                    g=gt.load_graph_from_csv(f'dataframes/{cl}/{cl}_{gn}_df_cluster0_{t}perc_tot_{flt}_linked.bed',skip_first=True,csv_options={'delimiter':'\t'})
                    g.save(f"graphs/{cl}/{cl}_{gn}_g_cluster0_{t}perc_tot_{flt}_linked.xml.gz")
                else:
                    print(cl,gn,t)
                    continue
        except FileNotFoundError:
            print(cl,gn)
            continue
