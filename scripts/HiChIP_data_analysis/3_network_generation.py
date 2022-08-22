import graph_tool.all as gt
import numpy as np
import pandas as pd


#1) obtain the list of bins of a specific chromatin iCD:
def neigh(df,param,anal,clique_tot):
    if clique_tot!=True:                  # If I am not interested in nodes of the whole "iCD" I return only the region
        return param                       #  where the promoter is located
    else:
        list_p=param                            #Else, for each region associated with the bin is stored
        analyzed=anal                         # in a list and one at a time are analyzed in order to find their neighbors
        todo=list(set(list_p)-set(analyzed))
        comparison=len(list_p)
        for r in todo:
            s_param=df[df['start']==r]
            e_param=df[df['end']==r]
            s_param_end=s_param['end'].to_list()
            e_param_start=e_param['start'].to_list()
            analyzed.append(r)
            for i in e_param_start:
                if i not in list_p:
                    list_p.append(i)
                else:
                    continue
            for i in s_param_end:
                if i not in list_p:
                    list_p.append(i)
                else:
                    continue
        if len(list_p)!=comparison:
            return neigh(df,list_p,analyzed,clique_tot)
        else:
            return list_p
            
#2) generation of a dataframe describing the iCD
def df_motif(df,list_of_regions):
    df_motifs=pd.DataFrame(columns=['start','end','contacts','mean_cov','peak1','peak2'])
    for r in list_of_regions:
        df_s=df[df['start']==r]
        df_e=df[df['end']==r]
        df_motifs=pd.concat([df_motifs,df_s,df_e])
    df_motifs=df_motifs[df_motifs.duplicated()==False]
    return df_motifs
  
#3) Order regions according to their order in genomic sequence  
def order_regions(df):
    reg=pd.DataFrame([x.split('-') for x in list(set(df['start'])|set(df['end']))])
    reg[1]=reg[1].astype(int)
    reg=reg.sort_values([0,1])
    reg[1]=reg[1].astype(str)
    reg.index=[x for x in range(len(reg.index))]
    reg['regions'] = reg[[0,1,2]].agg('-'.join, axis=1).astype('category')
    reg['codes']=[x for x in range(len(reg.index))]
    reg['chr']=reg[0].astype('category')
    reg=reg.iloc[:,-3:]
    return reg
 
#4) link sequential bins
def link_graph(df,reg):
    df['conn_type']=[0 for x in range(len(df.index))]
    for i in range(len(reg.index)):
        try:
            start_chr=reg['chr'].tolist()[i]
            end_chr=reg['chr'].tolist()[i+1]
            if start_chr == end_chr:
                region_df=pd.DataFrame(np.array([reg['regions'].tolist()[i],reg['regions'].tolist()[i+1],1,1,2,2,reg['codes'].tolist()[i],reg['codes'].tolist()[i+1],1]),index=['start','end','contacts','mean_cov','peak1','peak2','start_codes','end_codes','conn_type']).T
                df=pd.concat([df,region_df])
        except IndexError:
            continue
    df.index=[x for x in range(len(df.index))]
    df['start_codes']=df['start_codes'].astype(int)
    df['end_codes']=df['end_codes'].astype(int)
    df['contacts']=df['contacts'].astype(int)
    df['conn_type']=df['conn_type'].astype(int)
    df['start']=df['start'].tolist()
    df['start']=df['start'].astype('category')
    df['end']=df['end'].tolist()
    df['end']=df['end'].astype('category')
    return df


#5) Processing of HiChIP bed file
def prep_df(df,link,region,clique_tot):
    df['chr1']=df['chr1'].astype('str')
    df['s1']=df['s1'].astype('str')
    df['e1']=df['e1'].astype('str')
    df['chr2']=df['chr2'].astype('str')
    df['s2']=df['s2'].astype('str')
    df['e2']=df['e2'].astype('str')
    to_drop=['Bias1','Mapp1', 'GCContent1','RESites1','Bias2','Mapp2','GCContent2', 'RESites2','Dist','p','exp_cc_Bias', 'p_Bias','dbinom_Bias','P-Value_Bias','Q-Value_Bias']
    df=df.drop(to_drop,axis=1)
    df['start'] = df[['chr1','s1','e1']].agg('-'.join, axis=1).astype('category')
    df['end'] = df[['chr2','s2','e2']].agg('-'.join, axis=1).astype('category')
    df['contacts']=[i for i in df['cc'].tolist()]
    df['mean_cov']=(df['Coverage1']+df['Coverage2'])/2
    df['peak1']=df['isPeak1'].astype('category')
    df['peak2']=df['isPeak2'].astype('category')
    df=df.iloc[:,11:]    
    DF=pd.DataFrame(columns=['start','end','contacts','mean_cov','peak1','peak2'])
    if len(region)>0:
        for r in region:                                    # for each region all its connections are retrieved and iteratively
            list_of_neighbors=neigh(df,[r],[],clique_tot)     # stored in a df
            gf=df_motif(df,list_of_neighbors)
            DF=pd.concat([DF,gf])
        DF.index=[x for x in range(len(DF.index))]
        df=DF
    else:
        print('no regions')
    #create codes
    reg=order_regions(df)
    #<ssign the right code to each region
    start_codes=[]
    end_codes=[]
    for r in df['start']:
        if reg[reg['regions']==r].shape[0]==1:
            start_codes.append(reg[reg['regions']==r].iloc[0,1])
    for r in df['end']:
        if reg[reg['regions']==r].shape[0]==1:
            end_codes.append(reg[reg['regions']==r].iloc[0,1])    
    df['start_codes']=start_codes
    df['end_codes']=end_codes
    df['start']=df['start'].tolist()
    df['start']=df['start'].astype('category')
    df['end']=df['end'].tolist()
    df['end']=df['end'].astype('category')
    df['contacts']=df['contacts'].astype('int') #otherwise adj function does not work
    if link==True:
        df=link_graph(df,reg)
    df=df.drop_duplicates()
    df.index=[x for x in range(len(df.index))]
    return df,reg
    
#6) graph generation   
def get_graph(df,reg):
    g = gt.Graph(directed=False)
    g.add_vertex(reg.shape[0])
    if 'conn_type' in df.columns:            # if link between adjacent regions
        elist=[]
        for i in range(df.shape[0]):
            edge=(df.loc[i,'start_codes'],df.loc[i,'end_codes'],df.loc[i,'contacts'],df.loc[i,'conn_type'])
            elist.append(edge)
        eweight = g.new_ep("double")    #connections strength
        elayer = g.new_ep("int")     # layer (connection type: hichip=0, sequence=1)
        vtext=g.new_vp("int")     #vertex code
        g.add_edge_list(elist, eprops=[eweight, elayer])
        vtext.a=g.get_vertices()
        g.edge_properties['weight'] = eweight
        g.edge_properties['layer'] = elayer
        g.vertex_properties['id'] = vtext
    else:            # if link between adjacent regions
        elist=[]
        for i in range(df.shape[0]):
            edge=(df.loc[i,'start_codes'],df.loc[i,'end_codes'],df.loc[i,'contacts'])
            elist.append(edge)
        eweight = g.new_ep("double")    #connections strength
        vtext=g.new_vp("int")     #vertex code
        g.add_edge_list(elist, eprops=[eweight])
        vtext.a=g.get_vertices()
        g.edge_properties['weight'] = eweight
        g.vertex_properties['id'] = vtext
    return g
    
    
# Generation of MD6 total interactome network
cell='MD-6'
df=pd.read_csv(f'/home/HiChIP/{cell}-H3K27ac_10kb.interactions_FitHiC_Q0.01.bed', delimiter='\t')
df,reg=prep_df(df,False,[],True)
g_total=get_graph(df,reg)

# Generation of MD6 SOX9 iCD
cell='MD-6'
df=pd.read_csv(f'/home/HiChIP/{cell}-H3K27ac_10kb.interactions_FitHiC_Q0.01.bed', delimiter='\t')
region='chr17-70110000-70120000'                   # bin comprising SOX9 promoter
df,reg=prep_df(df,False,[region],True)
g_sox9=get_graph(df,reg)


