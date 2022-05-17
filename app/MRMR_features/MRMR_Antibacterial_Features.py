
import os
import sys
import pandas as pd
import re
from collections import Counter
import math
import numpy as np
import string
invalid_char = set(string.punctuation)


#########dpc ############
def ab_dpc_comp(df,q):
    result_DPC = []
    std = ["K","L","M","S"]
    std1 = ["Q","W","G","K","V"]
    zz = df.iloc[:,0]
    desired_DPC = ["KQ","KW","LW","MG","MK","SV"]
    header_list = []
    for i in desired_DPC:
        header_list.append(f"DPC_{i}")
    for i in range(0,len(zz)):
        result = []
        for j in std:
            for k in std1:
                count = 0
                temp = j+k
                if temp in desired_DPC:
                    for m3 in range(0,len(zz[i])-q):
                        b = zz[i][m3:m3+q+1:q]
                        if b == temp:
                            count += 1
                        composition = (count/(len(zz[i])-(q)))*100
                    result.append("%.2f" %composition)     
        result_DPC.append(result)  
        df_DPC = round(pd.DataFrame(result_DPC,columns = header_list),3) 
    return df_DPC

###################tpc ###############
def ab_tpc_comp(df):
    std = ["C","F","K","L","T","W"]
    std1 = ["A","F","K","L","M","Q","V","W"]
    std2 = ["A","D","I","K","L"]
    comb = ["CKI","FKK","KAA","KQL","KVL","KWK","LFK","LLD","TML","WKL"]
    header_list = []
    for i in comb:
        header_list.append(f"TPC_{i}")
    list_result = []
    zz = df.iloc[:,0]
    for i in range(0,len(zz)):
        list_2 = []
        for j in std:
            for k in std1:
                for m1 in std2:
                    count = 0
                    temp = j+k+m1
                    if temp in comb:
                        for m3 in range(0,len(zz[i])):
                            b = zz[i][m3:m3+3]
                            if b == temp:
                                count += 1
                            composition = (count/(len(zz[i])-2))*100
                        list_2.append("%.2f" %composition)
      
        list_result.append(list_2)
    df_TPC = round(pd.DataFrame(list_result,columns = header_list),3)
    return df_TPC



################## RRI #####################
def ab_RAAC(df):
    std = list("M")
    header_RRI = []
    for i in std:
        header_RRI.append(f"RRI_{i}")
    count = 0
    cc = []
    i = 0
    x = 0
    
    list_result = []
    for q in range(0,len(df)):
        list_1 = []
        while i < len(std):
            cc = []
            for j in df[0][q]:
                if j == std[i]:
                    count += 1
                    cc.append(count)
                else:
                    count = 0
            while x < len(cc) :
                if x+1 < len(cc) :
                    if cc[x]!=cc[x+1] :
                        if cc[x] < cc[x+1] :
                            cc[x]=0
                x += 1
            cc1 = [e for e in cc if e!= 0]
            cc = [e*e for e in cc if e != 0]
            zz= sum(cc)
            zz1 = sum(cc1)
            if zz1 != 0:
                zz2 = zz/zz1
            else:
                zz2 = 0
            
            list_1.append(zz2)
            i += 1
        i = 0
        list_result.append(list_1)
    df_RRI = round(pd.DataFrame(list_result,columns = header_RRI),3)
    return df_RRI


def ab_SE_residue_level(df):
    data = df[0].to_list()
    GH=[]
    for i in range(len(data)):
        my_list={'C':0}
        seq=data[i]
        num, length = Counter(seq), len(seq)
        num=dict(sorted(num.items()))
        C=list(num.keys())
        F=list(num.values())
        for key, value in my_list.items():
             for j in range(len(C)):
                if key == C[j]:
                    my_list[key] = round(((F[j]/length)* math.log(F[j]/length, 2)),3)
        GH.append(list(my_list.values()))
    header_SER = ["SER_C"]  
    df_SER = round(pd.DataFrame(GH,columns = header_SER),3)
    return df_SER


def ab_qos(df,gap,w=0.1):
    std = list("KW")
    mat1 = pd.read_csv(os.path.join(sys.path[0], "modal_csv", "Schneider-Wrede.csv"), index_col = 'Name')
    mat2 = pd.read_csv(os.path.join(sys.path[0],"modal_csv", "Grantham.csv"), index_col = 'Name')
    s1 = []
    s2 = []
    for i in range(0,len(df)):
        for n in range(1, gap+1):
            sum1 = 0
            sum2 = 0
            for j in range(0,(len(df[0][i])-n)):
                sum1 = sum1 + (mat1[df[0][i][j]][df[0][i][j+n]])**2
                sum2 = sum2 + (mat2[df[0][i][j]][df[0][i][j+n]])**2
            s1.append(sum1)
            s2.append(sum2)
    zz = pd.DataFrame(np.array(s1).reshape(len(df),gap))
    zz["sum"] = zz.sum(axis=1)
    zz2 = pd.DataFrame(np.array(s2).reshape(len(df),gap))
    zz2["sum"] = zz2.sum(axis=1)
    c1 = []
    c2 = []
    h1 = []
    h2 = []
    for aa in std:
        h1.append('QSO'+str(gap)+'_SC_' + aa)
    for aa in std:
        h2.append('QSO'+str(gap)+'_G_' + aa)
    for i in range(0,len(df)):
        AA = {}
        for j in std:
            AA[j] = df[0][i].count(j)
            c1.append(AA[j] / (1 + w * zz['sum'][i]))
            c2.append(AA[j] / (1 + w * zz2['sum'][i]))
    pp1 = np.array(c1).reshape(len(df),len(std))
    pp2 = np.array(c2).reshape(len(df),len(std))
    zz5 = round(pd.concat([pd.DataFrame(pp1, columns = h1),pd.DataFrame(pp2,columns = h2)],axis = 1),4)#pd.DataFrame(pp3, columns = h3),pd.DataFrame(pp4, columns = h4)], axis = 1),4)
    return zz5.loc[:,["QSO3_G_K","QSO3_SC_W"]]


##################CeTD##################
def ab_ctd(df):
    attr=pd.read_csv(os.path.join(sys.path[0],"modal_csv", "aa_attr_group.csv"), sep="\t")
    n = 0
    stt1 = []
    m = 1
    for i in range(0,len(attr)) :
        st =[]
        stt1.append([])
        for j in range(0,len(df)) :
            stt1[i].append([])
            for k in range(0,len(df[0][j])) :
                while m < 4 :
                    while n < len(attr.iloc[i,m]) :
                        if df[0][j][k] == attr.iloc[i,m][n] :
                            st.append(m)
                            stt1[i][j].append(m)
                        n += 2
                    n = 0
                    m += 1
                m = 1
#####################Composition##################### 
    std = [1,2,3]
    result_CTD = []
    for p in range (0,len(df)) :
        result = []
        for ii in range(0,len(stt1)) :
            for pp in std :
                count = 0
                for kk in stt1[ii][p] :
                    temp1 = kk
                    if temp1 == pp :
                        count += 1
                    composition = (count/len(stt1[ii][p]))*100
                result.append("%.2f"%composition)      
        result_CTD.append(result)    
        df_Comp = pd.DataFrame(result_CTD)    
    header1 = ['CeTD_HB','CeTD_VW','CeTD_PO','CeTD_PZ','CeTD_CH','CeTD_SS','CeTD_SA']
    head = []
    for i in header1:
        for j in range(1,4):
            head.append(i+str(j))
    df_Comp.columns = head

################ Transition #################
    tt = []
    tr=[]
    kk =0
    for ii in range(0,len(stt1)) :
        tt = []
        tr.append([])
        for p in range (0,len(df)) :
            tr[ii].append([])
            while kk < len(stt1[ii][p]) :
                if kk+1 <len(stt1[ii][p]):
                
                    tt.append(stt1[ii][p][kk])
                    tt.append(stt1[ii][p][kk+1])
                    tr[ii][p].append(stt1[ii][p][kk])
                    tr[ii][p].append(stt1[ii][p][kk+1])

                kk += 1
            kk = 0

    pp = 0
    xx = []
    xxx = []
    for mm in range(0,len(tr)) :
        xx = []
        xxx.append([])
        for nn in range(0,len(tr[mm])):
            xxx[mm].append([])
            while pp < len(tr[mm][nn]) :
                xx .append(tr[mm][nn][pp:pp+2])
                xxx[mm][nn].append(tr[mm][nn][pp:pp+2])
                pp+=2
            pp = 0
    std1 = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]]
    result_CTD = []
    for rr in range(0,len(df)) :
        result = []
        for qq in range(0,len(xxx)):
            for tt in std1 :
                count = 0
                for ss in xxx[qq][rr] :
                    temp2 = ss
                    if temp2 == tt :
                        count += 1
                result.append(count)         
        result_CTD.append(result)    
    df_Transition = pd.DataFrame(result_CTD)  
    head2 = []
    header2 = ['CeTD_11','CeTD_12','CeTD_1-3','CeTD_21','CeTD_22','CeTD_23','CeTD_31','CeTD_32','CeTD_33']
    for i in header2:
        for j in ('HB','VW','PO','PZ','CH','SS','SA'):
            head2.append(i+'_'+str(j))   
      
    df_Transition.columns = head2
    df_Comp_Transition = pd.concat([df_Comp,df_Transition],axis = 1)
    return df_Comp_Transition[["CeTD_23_HB","CeTD_31_PZ","CeTD_CH1","CeTD_SS1"]] 
    

############################# SEP ###############
import sys
PCP= pd.read_csv(os.path.join(sys.path[0], 'modal_csv', 'PhysicoChemical.csv'), header=None) 
headers = ['Positively charged','Negatively charged','Neutral charged','Polarity',
'Non polarity','Aliphaticity','Cyclic','Aromaticity','Acidicity','Basicity',
'Neutral (ph)','Hydrophobicity','Hydrophilicity','Neutral','Hydroxylic',
'Sulphur content','Secondary Structure(Helix)','Secondary Structure(Strands)','Secondary Structure(Coil)',
'Solvent Accessibilty (Buried)','Solvent Accesibilty(Exposed)','Solvent Accesibilty(Intermediate)','Tiny','Small','Large'];
def ab_encode(peptide):
    
    l=len(peptide);
    encoded=np.zeros(l);
    for i in range(l):
        if(peptide[i]=='A'):
            encoded[i] = 0;
        elif(peptide[i]=='C'):
            encoded[i] = 1;
        elif(peptide[i]=='D'):
            encoded[i] = 2;
        elif(peptide[i]=='E'):
            encoded[i] = 3;
        elif(peptide[i]=='F'):
            encoded[i] = 4;
        elif(peptide[i]=='G'):
            encoded[i] = 5;
        elif(peptide[i]=='H'):
            encoded[i] = 6;
        elif(peptide[i]=='I'):
            encoded[i] = 7;
        elif(peptide[i]=='K'):
            encoded[i] = 8;
        elif(peptide[i]=='L'):
            encoded[i] = 9;
        elif(peptide[i]=='M'):
            encoded[i] = 10;
        elif(peptide[i]=='N'):
            encoded[i] = 11;
        elif(peptide[i]=='P'):
            encoded[i] = 12;
        elif(peptide[i]=='Q'):
            encoded[i] = 13;
        elif(peptide[i]=='R'):
            encoded[i] = 14;
        elif(peptide[i]=='S'):
            encoded[i] = 15;
        elif(peptide[i]=='T'):
            encoded[i] = 16;
        elif(peptide[i]=='V'):
            encoded[i] = 17;
        elif(peptide[i]=='W'):
            encoded[i] = 18;
        elif(peptide[i]=='Y'):
            encoded[i] = 19;

    return encoded;

def ab_lookup(peptide,featureNum):
    l=len(peptide);
    peptide = list(peptide);
    out=np.zeros(l);
    peptide_num = ab_encode(peptide);
    
    for i in range(l):
        out[i] = PCP[peptide_num[i]][featureNum];
    return sum(out);

def ab_pcp(seq):
    l = len(seq);
    rows = PCP.shape[0]; 
    sequenceFeature = [];
    sequenceFeature.append(headers); 
    for i in range(l): 
        nfeatures = rows;
        sequenceFeatureTemp = [];
        for j in range(nfeatures): 
            featureVal = ab_lookup(seq[i],j)   
            if(len(seq[i])!=0):
                sequenceFeatureTemp.append(featureVal/len(seq[i]))
            else:
                sequenceFeatureTemp.append('NaN')
        sequenceFeature.append(sequenceFeatureTemp);
    return sequenceFeature;

def ab_phyChem(seq,mode='all',m=0,n=0):
    output = ab_pcp(seq);
    return output

def ab_sep_wp(df):
    seq=[]
    [seq.append(df.iloc[i][0]) for i in range(len(df))]
    comp = ab_phyChem(seq);
    new = [comp[i][0:25] for i in range(len(comp))]
    entropy  = [];
    for i in range(1,len(new)):
        seqEntropy = [];
        for j in range(len(new[i])):
            p = new[i][j]; 
            if((1-p) == 0. or p==0.):
                temp = 0;
            else:
                temp = -(p*math.log2(p)+(1-p)*math.log2(1-p));
            seqEntropy.append(round(temp,3));
        entropy.append(seqEntropy);
    Entropy_Dataframe = round(pd.DataFrame(entropy,columns = headers),3) 	
    
    desired_feature = ["Neutral"]
    return Entropy_Dataframe[desired_feature]    


