# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import pandas as pd
import numpy as np
#import requests
import pickle
import os
#from math import isnan
#from itertools import islice


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

def dataframe_uniprot(file_name,taxcatA, taxcatB): 
    """
    Read csv file and return pickled dataframe of UniprotID
    """
    df= pd.read_csv(file_name,chunksize=500)
    data = pd.DataFrame(columns=['idA','idB'])

    for df_chunk in df:
        df_uniprotID = df_pattern(df_chunk,taxcatA, taxcatB)  
        data = data.append(df_uniprotID, ignore_index=True)
    
    #drop duplicates
    data['check_string'] = data.apply(lambda row: ''.join(sorted([row['idA'], row['idB']])), axis=1)
    data.drop_duplicates('check_string')
    data.pop('check_string')
    
    file_name = "uniprot_ID_%s_%s.pkl" % (taxcatA, taxcatB)
    return data.to_pickle(file_name)

def dict_conversion_ID(filename):
    """
    Return a pickled dictionary for conversion purposes
    """
    key_value = (pd.read_csv(filename)).dropna()
    key_value_droppedna = key_value.dropna()
    
    #Seperate multiple entrez values seperated by ;
    key_value_droppedna['entrez']= key_value_droppedna['entrez'].apply(lambda x:
        list(filter(None, x.split(';'))))

    dict_conversion = (key_value_droppedna).set_index('uniprot').T.to_dict('list')
    pickle.dump( dict_conversion, open( "dict_conversion_ID.pkl", "wb" ) )

def df_uniprot_to_entrez(df,dict, taxcatA, taxcatB):
    
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

    file_name = "df_entrezID_%s_%s.pkl" % (taxcatA, taxcatB)
    return df_entrezID.to_pickle(file_name)

def make_dictionary_expression(directory):
    """
    Returns a dictionary that will enable determination of co-expression values
    """
    files = [i for i in os.listdir(directory)]
    
    conversion_dict = {} 
    for f in files:
        with open(os.path.join(directory,f)) as file_object:
            for lien in file_object:
             conversion_dict = {} 
             for lin in file_object:
                line = lin.split()
                key = (str(f),str(line[0]))
                if key not in conversion_dict:
                    conversion_dict[key] = str(line[2])
     
    pickle.dump(conversion_dict, open( "directorydict_expression.pkl", "wb" ) )
    
def make_dictionary_expression2(directory):
    """
    Returns a dictionary that will enable determination of co-expression values
    """
    files = [i for i in os.listdir(directory)]
    
    for f in files:
        with open(os.path.join(directory,f)) as file_object:
            name = directory+"/"+"dict_coexpression_" + f + ".pkl"
            rows = ( line.split('\t') for line in file_object ) 
            conversion_dict     = { (f,str(row[0])):str(row[2]) for row in rows}
            pickle.dump(conversion_dict, open( name, "wb" ) )


    
def entrez_to_coexpression(df, directory):
    """
    Returns pickled dataframe of coexpression values for
    a dataframe of entrez values via a coversion dictionary
    """
    df_coexpression = df.copy()
    df_coexpression['coexpression'] = np.empty((len(df), 0)).tolist()
    
    files = [i for i in os.listdir(directory)]
    
    for i, row in df_coexpression.iterrows():
        idA = (row['idA'])
        idB = (row['idB'])
        for elemA in idA:
            for elemB in idB:
                if elemA in files:
                    name = directory+"/"+"dict_coexpression_" + str(elemA) + ".pkl"
                    tup_tmp = (str(elemA), str(elemB))
                    dict = pickle.load( open( name, "rb" ) )
                    if tup_tmp in dict:
                        coexpression_value = float((dict.get(tup_tmp)).strip())
                        row['coexpression'].append(coexpression_value) 
        if (i%1000 == 0):
            print(i)
            df_coexpression.to_pickle("df_coexpression.pkl")
            
    #return df_coexpression       
            
#dataframe_uniprot("/media/sf_shared/ppi_data_2/ppi.csv", 'host', 'host')
#dict_conversion_ID("/media/sf_shared/ppi_data_2/uniprot2entrez.csv")
dict_conversion_ID = pickle.load( open( "dict_conversion_ID.pkl", "rb" ) )
#df_uniprotID= pd.read_pickle("uniprot_ID_host_host.pkl")
#df_uniprot_to_entrez(df_uniprotID, dict_conversion_ID, 'host', 'host')
#df_entrezID= pd.read_pickle("df_entrezID_host_host.pkl")
#df_entrezID.to_csv("/media/sf_shared/df_entrezID.csv")
    
#make_dictionary_expression2("/media/sf_shared/Co-expression_data/RNASEQo_expression_data_C")
#entrez_to_coexpression(df_entrezID,"/media/sf_shared/Co-expression_data/RNASEQ_Co_expression_data")
df_coexpression = pd.read_pickle("df_coexpression.pkl")
#df_coexpression.to_csv("/media/sf_shared/df_coexpression.csv")


#dataframe_uniprot("/media/sf_shared/ppi_data_2/ppi.csv", 'host', 'virus')
#df_uniprotID_host_virus= pd.read_pickle("uniprot_ID_host_virus.pkl")
#df_uniprotID_host_virus.to_csv("/media/sf_shared/df_uniprotID_host_virus.csv")
#df_uniprot_to_entrez(df_uniprotID_host_virus, dict_conversion_ID, 'host', 'virus')
df_entrezID_host_virus= pd.read_pickle("df_entrezID_host_virus.pkl")
df_entrezID_host_virus_reindexed = (df_entrezID_host_virus.set_index(['idB'])).sort_index()
#df_entrezID_host_virus.to_csv("/media/sf_shared/df_entrezID_host_virus.csv")

#dataframe_uniprot("/media/sf_shared/ppi_data_2/ppi.csv", 'host', 'bacteria')
#df_uniprotID_host_bacteria= pd.read_pickle("uniprot_ID_host_bacteria.pkl")
#df_uniprotID_host_bacteria.to_csv("/media/sf_shared/df_uniprotID_host_bacteria.csv")
#df_uniprot_to_entrez(df_uniprotID_host_bacteria, dict_conversion_ID, 'host', 'bacteria')
df_entrezID_host_bacteria= pd.read_pickle("df_entrezID_host_bacteria.pkl")
df_entrezID_host_bacteria_reindexed = (df_entrezID_host_bacteria.set_index(['idB'])).sort_index()
#df_entrezID_host_bacteria.to_csv("/media/sf_shared/df_entrezID_host_bacteria.csv")




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

