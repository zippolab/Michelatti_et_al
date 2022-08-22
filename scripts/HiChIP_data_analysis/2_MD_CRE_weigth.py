import pandas as pd
import graph_tool.all as gt
import random
from functions import network_generation  #folder containing functions for network generation

# THIS SCRIPT PRESENT THE ANALYSIS FOR MD CELLS TOTAL INTERACTIONS
# The same analysis have been performed for E-E, P-E, P-P interactions for tIMEC, XD and MD cells


# specific function for the analysis

def get_graph_2(df,reg):
    g = gt.Graph(directed=False)
    if 'conn_type' in df.columns:            # if link between adjacent regions
        elist=[]
        for i in range(df.shape[0]):
            edge=(df.loc[i,'start_codes'],df.loc[i,'end_codes'],df.loc[i,'contacts'],df.loc[i,'conn_type'])
            elist.append(edge)
        eweight = g.new_ep("double")    #connections strength
        elayer = g.new_ep("int")     # layer (connection type: hichip=0, sequence=1)
        vtext=g.new_vp("int")     #vertex code
        g.add_edge_list(elist, hashed=True,eprops=[eweight, elayer])
        vtext.a=g.get_vertices()
        g.edge_properties['weight'] = eweight
        g.edge_properties['layer'] = elayer
        g.vertex_properties['id'] = vtext
    else:            # if link between adjacent regions
        elist=pd.DataFrame([df['start_codes'].tolist(),df['end_codes'].tolist(),df['contacts'].tolist()]).T.to_numpy()
        eweight = g.new_ep("double")    #connections strength
        vtext=g.new_vp("int")     #vertex code
        g.add_edge_list(elist, hashed=True, eprops=[eweight])
        vtext.a=g.get_vertices()
        g.edge_properties['weight'] = eweight
        g.vertex_properties['id'] = vtext
    return g


# file loading
connections=pd.read_table('/home/HiChIP/MD-6-H3K27ac_10kb.interactions_FitHiC_Q0.01.bed')
connections['start']=connections['chr1']+'-'+connections['s1'].astype(str)+'-'+connections['e1'].astype(str)
connections['end']=connections['chr2']+'-'+connections['s2'].astype(str)+'-'+connections['e2'].astype(str)
connections['tot']=connections['start']+'_'+connections['end']

tot_intersection=pd.read_table('/home/HiChIP/references/HiChIP_MD_CRE_intersection.bed',header=None)

# input files generation
bin_MD=list(set(connections['start'])|set(connections['end']))
sample_size_MD=tot_intersection.shape[0]
print(sample_size_MD)                                # amount of bins that must be randomly sampled and the removed for each permutation

df,reg=network_generation.prep_df(connections,False,[],True)

# Average vertex degree w/o MD enriched CRE associated bins
to_remove=tot_intersection['bin'].tolist()
connections_filtered_i=df[~((df['start'].isin(to_remove))|(df['end'].isin(to_remove)))]
connections_filtered_i.index=[x for x in range(connections_filtered_i.shape[0])]
reg_filtered_i=reg[reg['regions'].isin(list(set(connections_filtered_i['start'])|set(connections_filtered_i['end'])))]
g_filtered_i_MD=get_graph_2(connections_filtered_i,reg_filtered_i)
print(gt.vertex_average(g_filtered_i_MD,deg='total')[0])          

# List of average vertex degree w/o randomly chosen bins (1000 permutations)
degree_dist_MD=[]
for i in range(1000):
    to_remove=sample(bin_MD,sample_size_MD)
    connections_filtered_i=df[~((df['start'].isin(to_remove))|(df['end'].isin(to_remove)))]
    connections_filtered_i.index=[x for x in range(connections_filtered_i.shape[0])]
    reg_filtered_i=reg[reg['regions'].isin(list(set(connections_filtered_i['start'])|set(connections_filtered_i['end'])))]
    g_filtered_i=get_graph_2(connections_filtered_i,reg_filtered_i)
    degree_dist_MD.append(gt.vertex_average(g_filtered_i,deg='total')[0])
    
    
pvalue=len([x for x in degree_dist_MD if x <gt.vertex_average(g_filtered_i_MD,deg='total')[0]])/len(degree_dist_MD)
