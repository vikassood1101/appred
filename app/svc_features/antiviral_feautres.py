import pandas as pd, os

################ ANTIVIRAL FEATURES #######################

############## bond ########

def AV_bond(df) :
    hy = []
    Du = []
    b2 = []
    b4 = []
    bb = pd.DataFrame()
    bonds=pd.read_csv(os.path.join("modal_csv", "bonds.csv"), sep = ",")
    for i in range(0,len(df)) :
        h = 0
        D = 0
        hy.append([i])
        Du.append([i])
        for j in range(0,len(df[0][i])) :
            temp = df[0][i][j]
            for k in range(0,len(bonds)) :
                if bonds.iloc[:,0][k] == temp :
                    h = h + bonds.iloc[:,2][k]
                    D = D + bonds.iloc[:,4][k]
        hy[i].append(h)
        Du[i].append(D)
    for m in range(0,len(df)) :
        b2.append(hy[m][1])
        b4.append(Du[m][1])
    bb["BTC_H"] = b2
    bb["BTC_D"] = b4
    return bb

###################### DDOR ##################
def AV_DDOR(df) :
    empty_list = []
    std = list("CFLRWY")
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


########################### CeTD ###############
def AV_ctd(df):
    attr=pd.read_csv(os.path.join("modal_csv", "aa_attr_group.csv"), sep="\t")
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
    features = ["CeTD_HB1",'CeTD_VW1','CeTD_PZ2','CeTD_CH1','CeTD_CH2','CeTD_CH3',"CeTD_SS2","CeTD_SA1","CeTD_SA2","CeTD_SA3"]
    df11 = df11.loc[:,features]
    return df11
    
################PAAC#####

def AV_paac_comp(df):
    std = list("GI")
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

################ APAAC##########
def AV_apaac_comp(df):
    std = list("GK")
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
