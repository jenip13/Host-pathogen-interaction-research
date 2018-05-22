# -*- coding: utf-8 -*-
"""
Functions to parse protein-protein interaction data

"""


import pandas as pd
import numpy as np
import pickle
import os


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
        
        if idA in files:
             name = directory+"/"+"dict_coexpression_" + str(idA) + ".pkl"
             tup_tmp = (str(idA), str(idB))
             dict1=  pickle.load( open( name, "rb" ) )
             if tup_tmp in dict1:
                 coexpression_value = float((dict1.get(tup_tmp)).strip())
                 row['coexpression'].append(coexpression_value) 
        
        if idB in files:
             name2 = directory+"/"+"dict_coexpression_" + str(idB) + ".pkl"
             tup_tmp2 = (str(idB), str(idA))
             dict2 = pickle.load( open( name2, "rb" ) )
             if tup_tmp2 in dict2:
                 coexpression_value2 = float((dict2.get(tup_tmp2)).strip())
                 if (coexpression_value != coexpression_value2 ):
                     row['coexpression'].append(coexpression_value2) 
    
    """
    
    for i, row in df_coexpression.iterrows():
        idA = (row['idA'])
        idB = (row['idB'])
        for elemA in list(idA):
            for elemB in list(idB):
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
       
    """     
    
    return df_coexpression       
            




