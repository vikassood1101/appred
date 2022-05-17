import pandas as pd, os
import numpy as np, sys
from .svc_features.amp_feautres import *
from .svc_features.anti_bacterial_feautres import*
from .svc_features.antiFungal_features import*
from .svc_features.antiviral_feautres import*
from .svc_features.nonAmp_features import *


################## SOC #########

def soc(df,gap):
    mat2 = pd.read_csv(os.path.join(sys.path[0], "modal_csv", "Grantham.csv"), index_col = 'Name')
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
    return pd.DataFrame(zz2)
    

########################## Function calling ################


def dataset_check(df,dataset):
    if dataset == "Antiviral":    
        df_AV_bond = AV_bond(df)
        df_AV_SOC = soc(df,3)
        df_AV_DDOR = AV_DDOR(df)  
        df_AV_PAAC = AV_paac_comp(df)
        df_AV_APAAC = AV_apaac_comp(df)
        df_AV_CTD = AV_ctd(df)
        df_result = pd.concat([df_AV_bond,df_AV_DDOR,df_AV_CTD,df_AV_PAAC,df_AV_APAAC,df_AV_SOC],axis = 1)
        
    
    elif dataset == "Antibacterial":
        df_AB_AAC = AB_aacp_comp(df)
        df_AB_DPC = AB_dpc_comp(df,1)
        df_AB_bond = AB_bond(df)
        df_AB_SOC = soc(df,3)
        df_AB_DDOR = AB_DDOR(df)  
        df_AB_PAAC = AB_paac_comp(df)
        df_AB_APAAC = AB_apaac_comp(df)
        df_AB_CTD = AB_ctd(df)
        df_result = pd.concat([df_AB_AAC,df_AB_DPC,df_AB_bond,df_AB_DDOR,df_AB_CTD,df_AB_PAAC,df_AB_APAAC,df_AB_SOC],axis = 1)
        
    
    elif dataset == "Antifungal":
        df_AF_Bond = AF_bond(df)
        df_AF_DDOR = AF_DDOR(df)
        df_AF_SOC = soc(df,3)
        df_AF_AAC = AF_aacp_comp(df)
        df_AF_DPC = AF_dpc_comp(df,1)
        df_AF_CTD = AF_ctd(df)
        df_AF_PAAC = AF_paac_comp(df)
        df_AF_APAAC = AF_apaac_comp(df)
        df_result = pd.concat([df_AF_AAC,df_AF_DPC,df_AF_Bond,df_AF_DDOR,df_AF_CTD,df_AF_PAAC,df_AF_APAAC,df_AF_SOC],axis = 1)
        
    
    elif dataset == "AMP":
        df_AMP_AAC = AMP_aacp_comp(df)
        df_AMP_bond = AMP_bond(df)
        df_AMP_SOC = soc(df,3)
        df_AMP_DDOR = AMP_DDOR(df)  
        df_AMP_APAAC = AMP_apaac_comp(df)
        df_AMP_CTD = AMP_ctd(df)
        df_result = pd.concat([df_AMP_AAC,df_AMP_bond,df_AMP_DDOR,df_AMP_CTD,df_AMP_APAAC,df_AMP_SOC],axis = 1)
       
    
    elif dataset == "NON_AMP":
        df_NON_AMP_aac = NON_AMP_aacp_comp(df)
        df_NON_AMP_bond = NON_AMP_bond(df)
        df_NON_AMP_SOC = soc(df,3)
        df_NON_AMP_DDOR = NON_AMP_DDOR(df)  
        df_NON_AMP_PAAC = NON_AMP_paac_comp(df)
        df_NON_AMP_APAAC = NON_AMP_apaac_comp(df)
        df_NON_AMP_CTD = NON_AMP_ctd(df)
        df_result = pd.concat([df_NON_AMP_aac,df_NON_AMP_bond,df_NON_AMP_DDOR,df_NON_AMP_CTD,df_NON_AMP_PAAC,df_NON_AMP_APAAC,df_NON_AMP_SOC],axis = 1)
    return df_result

def Predict_ALL_Features(dataset,method,sequence):
    
    
    if method == "predict":
        Fasta_df = [i[1].upper() for i in sequence ]
        Seq = [i[0] for i in sequence]
        df = pd.DataFrame(Fasta_df)
        df_result = dataset_check(df,dataset)
        
        
    elif method == "protein_scan":
        Fasta_df = [i[1].upper() for i in sequence ]
        Seq = [i[0] for i in sequence]
        df = pd.DataFrame(Fasta_df)
        df_result = dataset_check(df,dataset)
        
    
    elif method == "design": 
        Fasta_df = [i[1].upper() for i in sequence ]
        Seq = [i[0] for i in sequence]
        df = pd.DataFrame(Fasta_df)
        df_result = dataset_check(df,dataset)
    return Fasta_df,Seq,df_result
