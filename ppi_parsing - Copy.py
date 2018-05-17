# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import pandas as pd
import numpy as np
import requests
import mygene
import pickle
from itertools import islice


#Different interactions depending on type of species
def pattern(taxcatA, taxcatB):
    if taxcatA == 'virus' and taxcatB == 'virus':
        return 'virus_on_virus'
    if taxcatA == 'virus' and taxcatB == 'bacteria':
        return 'virus_on_bacteria'
    if taxcatA == 'bacteria' and taxcatB == 'bacteria':
        return 'bacteria_on_bacteria'
    if taxcatA == 'host' and taxcatB == 'host':
        return 'host_on_host'
    if taxcatA == 'host' and taxcatB == 'virus':
        return 'host_on_virus'
    if taxcatA == 'host' and taxcatB == 'bacteria':
        return 'host_on_bacteria'

# Return df containing list of interacting proteins betwen taxcatA and taxcatB
def df_pattern(df, taxcatA, taxcatB):
    # Create DataFrame where the interactions is only between taxcatA and taxcatB
    df_host_on_host = df[np.logical_and(df['taxcatA'] == taxcatA ,df['taxcatB'] == taxcatB)]
    df_unitprotID = df_host_on_host.loc[:,['idA','idB']]
    return df_unitprotID
    
def df_uniprot_to_entrez(df):
    
    mg = mygene.MyGeneInfo()
    list_uniprotID = list(zip(df['idA'],df['idB']))
    list_entrezID =[]
    for tup in list_uniprotID:
        try:
          if(mg.querymany(tup, scopes='uniprot', 
                                       species='human',entrezonly = True,
                                       returnall=True)).get('missing') == []:
              
              list_entrezID.extend(mg.querymany(tup, scopes='uniprot', 
                                       species='human',entrezonly = True,
                                       as_dataframe=True).loc[:,['_id']]['_id'].tolist())
        except:
               df_entrezID = pd.DataFrame(list_entrezID)
               df_entrezID.columns = ['idA', 'idB']
               
    
    df_entrezID = pd.DataFrame(list_entrezID)
    df_entrezID.columns = ['idA', 'idB']

    return df_entrezID

def coexpression(df):
    return

def process(file_name): 
    df= pd.read_csv("ppi.csv",chunksize=500)
    data = pd.DataFrame(columns=['idA','idB'])

    for df_chunk in df:
        df_uniprotID = df_pattern(df_chunk, 'host', 'host')  
        #df_entrezID = df_uniprot_to_entrez(df_uniprotID)
        #df_coexpress = coexpression(df_entrezID)
        data = data.append(df_uniprotID, ignore_index=True)
        #data.to_pickle("pickle.pkl")

    
    #drop duplicates
    data['check_string'] = data.apply(lambda row: ''.join(sorted([row['idA'], row['idB']])), axis=1)
    data.drop_duplicates('check_string')
    data.pop('check_string')
    return data

def df_uniprot_to_entrez2(df,dict):
    
    list_uniprotID = list(zip(df['idA'],df['idB']))
    list_entrezID =[]
    for row in list_uniprotID:
        if (row[0] in dict.keys() and row[1] in dict.keys()):
            sublist = dict.get(row[0])+dict.get(row[1])
            list_entrezID.append(sublist)
                         
    df_entrezID = pd.DataFrame(list_entrezID)
    df_entrezID.columns = ['idA', 'idB']

    return df_entrezID
            
#uniprot_ID = process("ppi.csv")
#(uniprot_ID).to_pickle("uniprot_ID.pkl")
#data = df_uniprot_to_entrez(uniprot_ID)
df_uniprotID= pd.read_pickle("uniprot_ID.pkl")
#df1 = data.iloc[:54000,:]
#df1.to_csv("uniprot_ID_part1.csv")
#df2 = data.iloc[54000:108000,:]
#df2.to_csv("uniprot_ID_part2.csv")
#df3 = data.iloc[108000:,:]
#df3.to_csv("uniprot_ID_part3.csv")

#key_value = pd.read_csv("uniprot_ID_merge.csv")
#dict_conversion = key_value.set_index('Uniprot').T.to_dict('list')
#pickle.dump( dict_conversion, open( "dict_conversion.pkl", "wb" ) )
#dict_conversion = pickle.load( open( "dict_conversion.pkl", "rb" ) )
#df_entrezID = df_uniprot_to_entrez2(df_uniprotID, dict_conversion)
#(df_entrezID ).to_pickle("df_entrezID.pkl")
df_entrezID= pd.read_pickle("df_entrezID.pkl")
df_entrezID.to_csv("df_entrezID.csv")

#df= pd.read_csv("ppi.csv",nrows=1000)
#df_uniprotID = df_pattern(df, 'host', 'host')  
#df_entrezID = df_uniprot_to_entrez(df_uniprotID)

with open('df_entrezID.csv', 'rb') as csvfile:
    lines_gen = islice(csvfile, 200)
    for line in lines_gen:
        print(type(line))
"""
<form action="/cgi-bin/coex_search.cgi" method="post">
		<input name="pair" value="1" type="hidden">
		<table>
		  <tbody><tr>
		    <th>
		      <img src="images/arrow_w.gif" class="middle" width="5" height="5">
		      <a href="javascript:;" onclick="setExample('box_EdgeAnnotation_pair','12462\t12465\n12461\t12465\n12461\t12468\n'); return false;" value="example"> example </a>
		    </th> 
		  </tr>
		  <tr>
		    <td>
		      <textarea name="genes" id="box_EdgeAnnotation_pair" cols="15" rows="9" onblur="displayQueryRule('box_EdgeAnnotation_pair')" onfocus="undiplayQueryRule('box_EdgeAnnotation_pair')" style="color: rgb(0, 0, 0);">		      </textarea>
		    </td>
		  </tr>
		</tbody></table>
		<div class="aCENTER">
		  <input value="submit" type="submit">
		</div>
	    </form>  
"""

payload = {
    'value': '111 11\n111 12',
    }

r = requests.post("http://coxpresdb.jp/cgi-bin/coex_search.cgi", 
                  data=payload)



