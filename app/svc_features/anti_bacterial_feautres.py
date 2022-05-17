
import pandas as pd, os, sys

########################### ANTIBACTERIAL FEATURES ################

#########################first feature  aac#####################

def AB_aacp_comp(df):
    std = list("K")
    AAC_Result = []
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
    return pd.DataFrame(AAC_Result)    

####################### second feature dpc ######################

def AB_dpc_comp(df,q):
    result_DPC = []
    std = ["L","R"]
    std1 = ["A","R"]
    zz = df.iloc[:,0]
    for i in range(0,len(zz)):
        result = []
        for j in std:
            for k in std1:
                count = 0
                temp = j+k
                if temp in ["LA","RR"]:
                    for m3 in range(0,len(zz[i])-q):
                        b = zz[i][m3:m3+q+1:q]
                        if b == temp:
                            count += 1
                        composition = (count/(len(zz[i])-(q)))*100
                    result.append("%.2f" %composition)     
        result_DPC.append(result)    
    return pd.DataFrame(result_DPC)    

############## third feature bond ######################

def AB_bond(df) :
    hy = []
    b2 = []
    bb = pd.DataFrame()
    bonds=pd.read_csv(os.path.join(sys.path[0], "modal_csv", "bonds.csv"), sep = ",")
    for i in range(0,len(df)) :
        h = 0
        hy.append([i])
        for j in range(0,len(df[0][i])) :
            temp = df[0][i][j]
            for k in range(0,len(bonds)) :
                if bonds.iloc[:,0][k] == temp :
                    h = h + bonds.iloc[:,2][k]
        hy[i].append(h)
        
    for m in range(0,len(df)) :
        b2.append(hy[m][1])
    bb["BTC_H"] = b2
    return bb
    

#################################  fourth feature  -----  ddor ---Distance Distribution  of residue   ----- ################

def AB_DDOR(df) :
    empty_list = []
    std = list("CDFLMRWY")
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
    df_DDOR = pd.DataFrame(empty_list)   
    return df_DDOR

########################### fifth feature CeTD ###############

def AB_ctd(df):
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
#####################Composition######################
    
    std = [1,2,3]
    result_CTD = []
    for p in range (0,len(df)) :
        result = []
        for ii in range(0,len(stt1)) :
            #for jj in stt1[ii][p]:
            for pp in std :
                count = 0
                for kk in stt1[ii][p] :
                    temp1 = kk
                    if temp1 == pp :
                        count += 1
                    composition = (count/len(stt1[ii][p]))*100
                result.append("%.2f"%composition)      
                
        result_CTD.append(result)    
        df11 = pd.DataFrame(result_CTD)    
 
    header1 = ['CeTD_HB','CeTD_VW','CeTD_PO','CeTD_PZ','CeTD_CH','CeTD_SS','CeTD_SA']
    head = []
    for i in header1:
        for j in range(1,4):
            head.append(i+str(j))
    df11.columns = head
    features = ["CeTD_VW3",'CeTD_PZ2','CeTD_CH3','CeTD_SS1','CeTD_SS3','CeTD_SA3']
    df11 = df11.loc[:,features]
    return df11
    

################## sixth feature PAAC3###################

def AB_paac_comp(df):
    std = list("K")
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
    return pd.DataFrame(PAAC_Result)    

####################### seventh feature APAAC3#############
def AB_apaac_comp(df):
    std = list("K")
    APAAC_Result = []
    
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
        APAAC_Result.append(Result)
    return pd.DataFrame(APAAC_Result)
