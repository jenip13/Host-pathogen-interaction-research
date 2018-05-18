# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import pandas as pd
import numpy as np
import requests
import pickle
import os
from math import isnan
from itertools import islice


def pattern(taxcatA, taxcatB):
    """
    Returns the type of interaction depending on the type of species
    """
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

def df_pattern(df, taxcatA, taxcatB):
    """
    Return df containing list of interacting proteins betwen taxcatA and taxcatB
    """
    # Create DataFrame where the interactions is only between taxcatA and taxcatB
    df_host_on_host = df[np.logical_and(df['taxcatA'] == taxcatA ,df['taxcatB'] == taxcatB)]
    df_unitprotID = df_host_on_host.loc[:,['idA','idB']]
    return df_unitprotID

def dataframe_uniprot(file_name): 
    """
    Read csv file and return pickled dataframe of UniprotID
    """
    df= pd.read_csv("ppi.csv",chunksize=500)
    data = pd.DataFrame(columns=['idA','idB'])

    for df_chunk in df:
        df_uniprotID = df_pattern(df_chunk, 'host', 'host')  
        data = data.append(df_uniprotID, ignore_index=True)
    
    #drop duplicates
    data['check_string'] = data.apply(lambda row: ''.join(sorted([row['idA'], row['idB']])), axis=1)
    data.drop_duplicates('check_string')
    data.pop('check_string')
    return data.to_pickle("uniprot_ID.pkl")

def dict_conversion(filename):
    """
    Return a pickled dictionary for conversion purposes
    """
    key_value = (pd.read_csv(filename)).dropna()
    key_value_droppedna = key_value.dropna()
    
    #Seperate multiple entrez values seperated by ;
    key_value_droppedna['entrez']= key_value_droppedna['entrez'].apply(lambda x:
        list(filter(None, x.split(';'))))

    dict_conversion = (key_value_droppedna).set_index('uniprot').T.to_dict('list')
    pickle.dump( dict_conversion, open( "dict_conversion.pkl", "wb" ) )

def df_uniprot_to_entrez(df,dict):
    
    """
    Convert Dataframe with UniprotID to a pickled Dataframe with EntrezID
    via a conversion dictionary
    """
    
    list_uniprotID = list(zip(df['idA'],df['idB']))
    list_entrezID =[]
    for row in list_uniprotID:
        if (row[0] in dict.keys() and row[1] in dict.keys()):
            idA = dict.get(row[0])
            idB = dict.get(row[1])
            """
            if len(idB) >0:
                for element in idB:
                    sublist = idA + [element]
                    list_entrezID.append(sublist)    
            else:
            """
            sublist = idA + idB
            list_entrezID.append(sublist)
                         
    df_entrezID = pd.DataFrame(list_entrezID)
    df_entrezID.columns = ['idA', 'idB']

    return df_entrezID.to_pickle("df_entrezID.pkl")

def make_dictionary_expression(filepath):
    """
    Returns a dictionary
    """
    
def entrez_to_coexpression(df, dict):
    """
    Returns pickled dataframe of coexpression values for
    a dataframe of entrez values via a coversion dictionary
    """
    directory = "D:\shared\Co-expression_data\RNASEQ_Co_expression_data" 
    file_extension = '.csv' 
    
    conversion_dict = {} 

    files = [i for i in os.listdir(directory) if os.path.splitext(i)[1] == file_extension]

# Iterate over your txt files
for f in txt_files:
    # Open them and assign them to file_dict
    with open(os.path.join(direc,f)) as file_object:
        file_dict[f] = file_object.read()


    for row in df.itterows():
         for :
             for :
                 
            
#dataframe_uniprot("ppi.csv")
#dict_conversion("uniprot2entrez.csv")

dict_conversion = pickle.load( open( "dict_conversion.pkl", "rb" ) )
df_uniprotID= pd.read_pickle("uniprot_ID.pkl")

#df_uniprot_to_entrez(df_uniprotID, dict_conversion)
df_entrezID= pd.read_pickle("df_entrezID.pkl")
#df_entrezID.to_csv("df_entrezID.csv")




"""
with open('df_entrezID.csv', 'rb') as csvfile:
    lines_gen = islice(csvfile, 200)
    for line in lines_gen:
        print(type(line))

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


payload = {
    'value': '111 11\n111 12',
    }

r = requests.post("http://coxpresdb.jp/cgi-bin/coex_search.cgi", 
                  data=payload)

"""

