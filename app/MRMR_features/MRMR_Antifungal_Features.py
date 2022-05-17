import pandas as pd
from collections import Counter
import math, os, sys
import numpy as np
import string
invalid_char = set(string.punctuation)

######aac ##########
def af_aacp_comp(df):
    std = list("C")
    AAC_Result = []
    Header_name = []
    for i in std:
        Header_name.append(f"AAC_{i}")
    zz = df.iloc[:,0]
    for j in zz:
        Result = []
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            Result.append(composition) 
        AAC_Result.append(Result)
    df_AAC = round(pd.DataFrame(AAC_Result,columns = Header_name),3) 
    return df_AAC

#########dpc ############
def af_dpc_comp(df,q):
    result_DPC = []
    std = ["A","F","K","L","N","S","W"]
    std1 = ["C","K","P","R","V","Y"]
    zz = df.iloc[:,0]
    desired_DPC = ["AR","FP","KC","KK","LC","NY","SV","WK"]
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

############### tpc #################
def af_tpc_comp(df):
    std = ["C","E","F","L","V"]
    std1 = ["K","F"]
    std2 = ["I","V","K","P"]
    comb = ["CKI","EKV","FKK","LKV","VFP"]
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
def af_DDOR(df) :
    empty_list = []
    std = list("CKN")
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

################## RRI #####################
def af_RAAC(df):
    std = list("C")
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


############### QOS ################
def af_qos(df,gap,w=0.1):
    std = list("KC")
    mat1 = pd.read_csv(os.path.join(sys.path[0],"modal_csv", "Schneider-Wrede.csv"), index_col = 'Name')
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
    return zz5.loc[:,["QSO3_G_K","QSO3_SC_C"]]

###### SER ###########

def af_SE_residue_level(df):
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

################ APAAAC ############

def af_apaac_1(df,lambdaval,w=0.05):
    data1 = pd.read_csv(os.path.join(sys.path[0],"modal_csv","data"), sep = "\t")
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
    pp = pd.DataFrame()
    ee = []
    for k in range(0,len(df)):
        cc = [] 
        for n in range(1,lambdaval+1):
            for b in range(0,len(zz)):
                cc.append(sum([zz.loc[b][aa[df[0][k][p]]] * zz.loc[b][aa[df[0][k][p + n]]] for p in range(len(df[0][k]) - n)]) / (len(df[0][k]) - n))
                qq = pd.DataFrame(cc)
        pseudo = [(w * p) / (1 + w * sum(cc)) for p in cc]
        ee.append(pseudo)
        ii = round(pd.DataFrame(ee, columns = head),4)
    return pd.DataFrame(ii["APAAC3_HB_lam2"],columns=["APAAC3_HB_lam2"] )     

######################## SOC ##################
def af_soc(df,gap):
    mat2 = pd.read_csv(os.path.join(sys.path[0],"modal_csv", "Grantham.csv"), index_col = 'Name')
    h2 = []
    for n in range(1, gap + 1):
        h2.append('G' + str(n))
    h2 = ['SOC'+str(gap)+'_'+sam for sam in h2]
    s2 = []
    for i in range(0,len(df)):
        for n in range(1, gap+1):
            sum2 =0
            sum3 =0
            for j in range(0,(len(df[0][i])-n)):
                sum2 = sum2 + (mat2[df[0][i][j]][df[0][i][j+n]])**2
                sum3 = sum2/(len(df[0][i])-n)
            s2.append(sum3)
    zz2 = np.array(s2).reshape(len(df),gap)
    return pd.DataFrame(zz2,columns= h2)["SOC3_G1"]
    
##################### PAAC #####################
def af_paac_comp(df):
    std = list("C")
    header = []
    for i in std:
        header.append(f"PAAC3_{i}")
    PAAC_Result = []
    zz = df.iloc[:,0]
    for j in zz:
        Result = []
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            Result.append(composition) 
        PAAC_Result.append(Result)
    return pd.DataFrame(PAAC_Result,columns= header)        

######################### ATC ###################
def af_atc(df):
    atom=pd.read_csv(os.path.join(sys.path[0],"modal_csv", "atom.csv"),header=None)
    at=pd.DataFrame()
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
    return final[["ATC_H"]]

