import pandas as pd
import math, os, sys
import numpy as np
import string
invalid_char = set(string.punctuation)

################ dpc ######################
def amp_dpc_comp(df,q):
    result_DPC = []
    std = ["L","S"]
    std1 = ["K","V"]
    zz = df.iloc[:,0]
    desired_DPC = ["LK","SV"]
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


######### TPC #######
def amp_tpc_comp(df):
    std = ["A","G"]
    std1 = ["G","L","V","K"]
    std2 = ["A","L","W","K"]
    comb = ["AGK","ALW","AVL","GKA"]
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

############## DDOR #######################

def amp_DDOR(df) :
    empty_list = []
    std = list("DEL")
    header_DDOR = []
    for i in std:
        header_DDOR.append(f"DDOR_{i}")
    for i in range(0,len(df)):
        list_1 = []
        s = df[0][i]
        p = s[::-1]
        for j in std:
            zz = ([pos for pos, char in enumerate(s) if char == j])
            pp = ([pos for pos, char in enumerate(p) if char == j])
            ss = []
            for i in range(0,(len(zz)-1)):
                ss.append(zz[i+1] - zz[i]-1)
            if zz == []:
                ss = []
            else:
                ss.insert(0,zz[0])
                ss.insert(len(ss),pp[0])
            cc1=  (sum([e for e in ss])+1)
            cc = sum([e*e for e in ss])
            zz2 = cc/cc1
            list_1.append("%.2f"%zz2)
        empty_list.append(list_1)
    df_DDOR = round(pd.DataFrame(empty_list,columns = header_DDOR),3)
    return df_DDOR

################# QSO ############
def amp_qos(df,gap,w=0.1):
    std = list("KL")
    mat1 = pd.read_csv(os.path.join(sys.path[0], "modal_csv", "Schneider-Wrede.csv"), index_col = 'Name')
    mat2 = pd.read_csv(os.path.join(sys.path[0], "modal_csv", "Grantham.csv"), index_col = 'Name')
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
    return zz5.loc[:,["QSO3_G_K","QSO3_G_L","QSO3_SC_K"]]

############## PAAAC ###########

def amp_val(AA_1, AA_2, aa, mat):
    return sum([(mat[i][aa[AA_1]] - mat[i][aa[AA_2]]) ** 2 for i in range(len(mat))]) / len(mat)
def amp_paac_1(df,lambdaval,w=0.05):
    std_P = list("ACDEFGHIKLMNPQRSTVWY")
    data1 = pd.read_csv(os.path.join(sys.path[0], "modal_csv","data"), sep = "\t")
    dd = []
    cc = []
    pseudo = []
    aa = {}
    for i in range(len(std_P)):
        aa[std_P[i]] = i
    for i in range(0,3):
        mean = sum(data1.iloc[i][1:])/20
        rr = math.sqrt(sum([(p-mean)**2 for p in data1.iloc[i][1:]])/20)
        dd.append([(p-mean)/rr for p in data1.iloc[i][1:]])
    head = []
    for n in range(1, lambdaval + 1):
        head.append('_lam' + str(n))
    head = ['PAAC'+str(lambdaval)+sam for sam in head]
    ee = []
    for k in range(0,len(df)):
        cc = []
        pseudo1 = [] 
        for n in range(1,lambdaval+1):
            cc.append(sum([amp_val(df[0][k][p], df[0][k][p + n], aa, dd) for p in range(len(df[0][k]) - n)]) / (len(df[0][k]) - n))
        pseudo = pseudo1 + [(w * p) / (1 + w * sum(cc)) for p in cc]
        ee.append(pseudo)
        ii = round(pd.DataFrame(ee, columns = head),4)
    return pd.DataFrame(ii["PAAC3_lam2"],columns=["PAAC3_lam2"] ) 

#################### APAAC #####################

def amp_apaac_1(df,lambdaval,w=0.05):
    data1 = pd.read_csv(os.path.join(sys.path[0], "modal_csv", "data"), sep = "\t")
    std_AP = list("ACDEFGHIKLMNPQRSTVWY")
    dd = []
    cc = []
    pseudo = []
    aa = {}
    for i in range(len(std_AP)):
        aa[std_AP[i]] = i
    for i in range(0,3):
        mean = sum(data1.iloc[i][1:])/20
        rr = math.sqrt(sum([(p-mean)**2 for p in data1.iloc[i][1:]])/20)
        dd.append([(p-mean)/rr for p in data1.iloc[i][1:]])
        zz = pd.DataFrame(dd)
    head = []
    for n in range(1, lambdaval + 1):
        for e in ('HB','HL','SC'):
            head.append(e+'_lam' + str(n))
    head = ['APAAC'+str(lambdaval)+'_'+sam for sam in head]
    
    ee = []
    for k in range(0,len(df)):
        cc = [] 
        for n in range(1,lambdaval+1):
            for b in range(0,len(zz)):
                cc.append(sum([zz.loc[b][aa[df[0][k][p]]] * zz.loc[b][aa[df[0][k][p + n]]] for p in range(len(df[0][k]) - n)]) / (len(df[0][k]) - n))
                
        pseudo = [(w * p) / (1 + w * sum(cc)) for p in cc]
        ee.append(pseudo)
        ii = round(pd.DataFrame(ee, columns = head),4)
    return pd.DataFrame(ii["APAAC3_HL_lam2"],columns=["APAAC3_HL_lam2"] )     


################### ATC ##########################    
def amp_atc(df):
    atom=pd.read_csv(os.path.join(sys.path[0], "modal_csv","atom.csv"),header=None)
    i = 0
    C_atom = []
    H_atom = []
    N_atom = []
    O_atom = []
    S_atom = []
    while i < len(atom):
        C_atom.append(atom.iloc[i,1].count("C"))
        H_atom.append(atom.iloc[i,1].count("H"))
        N_atom.append(atom.iloc[i,1].count("N"))
        O_atom.append(atom.iloc[i,1].count("O"))
        S_atom.append(atom.iloc[i,1].count("S"))
        i += 1
    atom["C_atom"]=C_atom
    atom["O_atom"]=O_atom
    atom["H_atom"]=H_atom
    atom["N_atom"]=N_atom
    atom["S_atom"]=S_atom
    
    count_C = 0
    count_H = 0
    count_N = 0
    count_O = 0
    count_S = 0
    
    i1 = 0
    j = 0
    k = 0
    C_ct = []
    H_ct = []
    N_ct = []
    O_ct = []
    S_ct = []
    while i1 < len(df) :
        while j < len(df[0][i1]) :
            while k < len(atom) :
                if df.iloc[i1,0][j]==atom.iloc[k,0].replace(" ","") :
                    count_C = count_C + atom.iloc[k,2]
                    count_H = count_H + atom.iloc[k,3]
                    count_N = count_N + atom.iloc[k,4]
                    count_O = count_O + atom.iloc[k,5]
                    count_S = count_S + atom.iloc[k,6]
                
                k += 1
            k = 0
            j += 1
        C_ct.append(count_C)
        H_ct.append(count_H)
        N_ct.append(count_N)
        O_ct.append(count_O)
        S_ct.append(count_S)
        count_C = 0
        count_H = 0
        count_N = 0
        count_O = 0
        
        j = 0
        i1 += 1
    df["C_count"]=C_ct
    df["H_count"]=H_ct
    df["N_count"]=N_ct
    df["O_count"]=O_ct
    df["S_count"]=S_ct

    ct_total = []
    m = 0
    while m < len(df) :
        ct_total.append(df.iloc[m,1] + df.iloc[m,2] + df.iloc[m,3] + df.iloc[m,4] + df.iloc[m,5])
        m += 1
    df["count"]=ct_total

    final = pd.DataFrame()
    n = 0
    
    C_p = []
    H_p = []
    N_p = []
    O_p = []
    S_p = []
    while n < len(df):
        C_p.append((df.iloc[n,1]/df.iloc[n,6])*100)
        H_p.append((df.iloc[n,2]/df.iloc[n,6])*100)
        N_p.append((df.iloc[n,3]/df.iloc[n,6])*100)
        O_p.append((df.iloc[n,4]/df.iloc[n,6])*100)
        S_p.append((df.iloc[n,5]/df.iloc[n,6])*100)
        n += 1
    final["ATC_C"] = C_p
    final["ATC_H"] = H_p
    final["ATC_N"] = N_p
    final["ATC_O"] = O_p
    final["ATC_S"] = S_p
    return final[["ATC_H","ATC_N"]]

############### SEP ################
import sys
PCP= pd.read_csv(os.path.join(sys.path[0],'modal_csv', 'PhysicoChemical.csv'), header=None) 
headers = ['Positively charged','Negatively charged','Neutral charged','Polarity',
'Non polarity','Aliphaticity','Cyclic','Aromaticity','Acidicity','Basicity',
'Neutral (ph)','Hydrophobicity','Hydrophilicity','Neutral','Hydroxylic',
'Sulphur content','Secondary Structure(Helix)','Secondary Structure(Strands)','Secondary Structure(Coil)',
'Solvent Accessibilty (Buried)','Solvent Accesibilty(Exposed)','Solvent Accesibilty(Intermediate)','Tiny','Small','Large'];
def amp_encode(peptide):
    
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

def amp_lookup(peptide,featureNum):
    l=len(peptide);
    peptide = list(peptide);
    out=np.zeros(l);
    peptide_num = amp_encode(peptide);
    
    for i in range(l):
        out[i] = PCP[peptide_num[i]][featureNum];
    return sum(out);

def amp_pcp(seq):
    l = len(seq);
    rows = PCP.shape[0]; 
    sequenceFeature = [];
    sequenceFeature.append(headers); 
    for i in range(l): 
        nfeatures = rows;
        sequenceFeatureTemp = [];
        for j in range(nfeatures): 
            featureVal = amp_lookup(seq[i],j)   
            if(len(seq[i])!=0):
                sequenceFeatureTemp.append(featureVal/len(seq[i]))
            else:
                sequenceFeatureTemp.append('NaN')
        sequenceFeature.append(sequenceFeatureTemp);
    return sequenceFeature;

def amp_phyChem(seq,mode='all',m=0,n=0):
    output = amp_pcp(seq);
    return output

def amp_sep_wp(df):
    seq=[]
    [seq.append(df.iloc[i][0]) for i in range(len(df))]
    comp = amp_phyChem(seq);
    new = [comp[i][0:25] for i in range(len(comp))]
    entropy  = [];
    #entropy.append(headers[0:25])
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
    
    desired_feature = ["Negatively charged","Neutral","Positively charged","Polarity","Hydroxylic"]
    return Entropy_Dataframe[desired_feature]    

################# ctd #####################

def amp_ctd(df):
    attr=pd.read_csv(os.path.join(sys.path[0], "modal_csv", "aa_attr_group.csv"), sep="\t")
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
   
################# Distribution ################    
    
    c_22 = []
    for x in range(0,len(stt1)) :
        c_22.append([])
        k = 0
        j = 0
        for y in range(0,len(stt1[x])):
            
            c_22[x].append([])
            for i in range(1,4) :
                cc = []
                c1 = [index for index,value in enumerate(stt1[x][y]) if value == i]
                c_22[x][y].append(c1)
    cc = []

    for ss in range(0,len(df)):
        for uu in range(0,len(c_22)):
            for mm in range(0,3):
                for ee in range(0,101,25):
                    k = (ee*(len(c_22[uu][ss][mm])))/100
                    cc.append(math.floor(k))
                    
    x = [cc[i:i+105] for i in range(0, len(cc), 105)] 
    df_Distribution = pd.DataFrame(x)
    
    head3 = []
    header3 = ['CeTD_0_p','CeTD_25_p','CeTD_50_p','CeTD_75_p','CeTD_100_p']
    header4 = ['HB','VW','PO','PZ','CH','SS','SA']
    for j in range(1,4):
        for k in header4:
            for i in header3:
                head3.append(i+'_'+k+str(j))
     
    df_Distribution.columns = head3     
    df_Comp_distribution = pd.concat([df_Comp,df_Distribution],axis = 1)
    
    return df_Comp_distribution[["CeTD_25_p_SA3","CeTD_50_p_SA3","CeTD_SA1","CeTD_SA3"]]

