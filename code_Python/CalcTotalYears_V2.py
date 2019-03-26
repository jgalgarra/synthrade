# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 12:27:54 2017

@author: maryluz.mouronte (email: maryluz.mouronte@ufv.es)
"""

"""
Editor de Spyder

Este es un archivo temporal
"""

import pandas as pd
import numpy as np
import os
import os.path as path
import numpy as np
import pandas as pd
import networkx as nx
import operator
import random
from random import randint
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from os import listdir
import shutil
import time
from collections import OrderedDict
import copy
import glob
import dask.dataframe as dd



def LeerFicheros():
 
    
    Header = ['year','origin', 'dest', 'hs07', 'export_val','import_val']
    #fileName = '..\datos\year_origin_destination_hs07_6.csv'
    #fileNameIn = '..\datos\year_origin_destination_sitc_rev2.tsv'
    fileNameIn = '..\datos\year_origin_destination_hs07_6.tsv'
    #fileNameIn = '..\datos\part1.tsv'
    pathFileIn = os.getcwd() + "\\" + fileNameIn

    fileNameInP = '..\datos\Pyear_origin_destination_hs07_6.tsv'
    #fileNameInP = '..\datos\Ppart1.tsv'
    pathFileInP = os.getcwd() + "\\" + fileNameInP
    if path.exists(pathFileInP):
        os.remove(pathFileInP)  


    print('Leyendo fichero') 
    
    
    r = open(pathFileIn,"r")
    w = open (pathFileInP, "w")
    l = r.readline()
    cont = 1
    while l:
      year = l.split('\t')[0]
     
      
      
      if (year.isdigit()):
          if (int(year) in Years):
              l=l.replace('\t',';')
              w.write(l) 
          
      elif (not(year.isdigit())):
          l=l.replace('\t',';')
          w.write(l) #Escribe la cabecera
      print(cont)  
      
      l= r.readline()
      cont=cont+1
          
          
    r.close()
    w.close()

    datosRed = pd.read_table(pathFileInP, 'engine=python', delimiter=';', header=0, encoding = "ISO-8859-1", names=Header)
     
    #print(datosRed.head(10))
    print("Esto es")
    print(datosRed)
    for i in Years:
        print("Buildind year")
        print(i)
        datosRedRes = datosRed[(datosRed.year == i)]
        datosRedRes = datosRedRes.drop(labels="hs07", axis=1)
        print(datosRedRes.head(10))
        #datosRedRes = datosRedRes.groupby(by=['year','origin', 'dest'])['export_val','import_val'].sum()
        datosRedRes=datosRedRes.fillna(0)
        print(datosRedRes.head(10))
        print(datosRedRes.dtypes)
        datosRedRes = datosRedRes.groupby(['year', 'origin', 'dest'])['export_val','import_val'].sum()
        #datosRedRes = datosRedRes.groupby(['year', 'origin', 'dest'])[['export_val','import_val']].sum()
        print("Despues de agrupar")
        print(datosRedRes.head(10))
        fileName = '..\datos\TodosYR'+str(i)+'.csv'
        pathFile = os.getcwd() + "\\" + fileName
        Header = ['year', 'origin', 'dest', 'export_val', 'import_val']  
        
       
        datosRedRes.reset_index().to_csv(pathFile, index=None, mode='a', sep=';', columns=Header)
     


Years = [2015, 2016, 2017]   
print("Building file")
LeerFicheros()  








