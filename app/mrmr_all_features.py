from .MRMR_features.MRMR_AMP_Features import*
from .MRMR_features.MRMR_Antibacterial_Features import*
from .MRMR_features.MRMR_Antifungal_Features import*
from .MRMR_features.MRMR_Antiviral_Features import*
from .MRMR_features.MRMR_NonAMP_Features import*

def mrmr_dataset_check(df,dataset):
    if dataset == "Antiviral":    
        df_AAC = av_aacp_comp(df)
        df_DPC = av_dpc_comp(df,1)
        df_TPC = av_tpc_comp(df)
        df_ATC = av_atc(df)
        df_DDOR = av_DDOR(df)
        df_RRI = av_RAAC(df)
        df_SER = av_SE_residue_level(df)
        df_SEP= av_sep_wp(df)
        df_CTD = av_ctd(df)
        df_PAAC = av_paac_1(df,3)
        df_APAAC = av_apaac_1(df,3)
        df_QSO = av_qos(df,3)
        df_SOC = av_soc(df,3)
        df_result = pd.concat([df_AAC,df_DPC,df_TPC,df_ATC,df_DDOR,df_RRI,df_SER,df_SEP,df_CTD,df_PAAC,df_APAAC,df_QSO,df_SOC],axis = 1)
        column_names = ["ATC_H","RRI_K","DDOR_D","APAAC3_HL_lam2","CeTD_23_CH","ATC_N","DPC_WK","CeTD_CH3","PAAC3_lam2","SER_K","DPC_LK","DPC_GK","Neutral",
        "DDOR_E","AAC_G","QSO3_G_K","CeTD_SA1","SOC3_SC2","Positively charged","Negatively charged","TPC_FKK","CeTD_23_PO","AAC_K","Polarity","QSO3_SC_D"]
        df_result = df_result.reindex(columns=column_names)
    
    elif dataset == "Antibacterial":
        df_DPC = ab_dpc_comp(df,1)
        df_TPC = ab_tpc_comp(df)
        df_RRI = ab_RAAC(df)
        df_SER = ab_SE_residue_level(df)
        df_SEP= ab_sep_wp(df)
        df_CTD = ab_ctd(df)
        df_QSO = ab_qos(df,3)
        
        df_result = pd.concat([df_QSO,df_RRI,df_DPC,df_TPC,df_SER,df_CTD,df_SEP],axis = 1)
        column_names = ['CeTD_23_HB','CeTD_SS1','RRI_M','DPC_KQ', 'DPC_KW','TPC_KAA','QSO3_G_K','DPC_SV','TPC_LFK','TPC_CKI','Neutral','CeTD_31_PZ',
        'TPC_KWK', 'TPC_KVL','SER_C','TPC_KQL','CeTD_CH1','TPC_FKK','DPC_LW','DPC_MK','QSO3_SC_W','DPC_MG','TPC_TML','TPC_WKL','TPC_LLD']
        df_result = df_result.reindex(columns=column_names)
        
    
    elif dataset == "Antifungal":
        df_AAC = af_aacp_comp(df)
        df_DPC = af_dpc_comp(df,1)
        df_TPC = af_tpc_comp(df)
        df_ATC = af_atc(df)
        df_DDOR = af_DDOR(df)
        df_RRI = af_RAAC(df)
        df_SER = af_SE_residue_level(df)
        df_PAAC = af_paac_comp(df)
        df_APAAC = af_apaac_1(df,3)
        df_QSO = af_qos(df,3)
        df_SOC = af_soc(df,3)
        
        df_result = pd.concat([df_ATC,df_PAAC,df_AAC,df_DPC,df_TPC,df_DDOR,df_RRI,df_SER,df_QSO,df_APAAC,df_SOC],axis = 1)
        column_names = ['QSO3_SC_C','APAAC3_HB_lam2','DPC_NY','DPC_WK','RRI_C','DPC_SV','DDOR_K','DPC_KC','TPC_CKI','SER_C','DPC_KK','DPC_FP',
        'SOC3_G1','DDOR_N','DPC_LC','DDOR_C','DPC_AR','TPC_LKV','QSO3_G_K','AAC_C','TPC_VFP','TPC_FKK','ATC_H','TPC_EKV','PAAC3_C']
        df_result = df_result.reindex(columns=column_names)
        
    
    elif dataset == "AMP":
        df_DPC = amp_dpc_comp(df,1)
        df_TPC = amp_tpc_comp(df)
        df_ATC = amp_atc(df)
        df_DDOR = amp_DDOR(df)
        df_SEP= amp_sep_wp(df)
        df_CTD = amp_ctd(df)
        df_PAAC = amp_paac_1(df,3)
        df_APAAC = amp_apaac_1(df,3)
        df_QSO = amp_qos(df,3)
        df_result = pd.concat([df_CTD,df_ATC,df_DPC,df_DDOR,df_TPC,df_QSO,df_PAAC,df_APAAC,df_SEP],axis = 1)
        column_name = ['ATC_N','Hydroxylic','CeTD_SA1','QSO3_G_K','Neutral','DDOR_L','TPC_AGK','CeTD_50_p_SA3','PAAC3_lam2',
        'Negatively charged','QSO3_SC_K','DPC_SV','TPC_AVL','TPC_ALW','Polarity','DPC_LK','Positively charged',
        'CeTD_25_p_SA3','QSO3_G_L','APAAC3_HL_lam2','ATC_H','DDOR_E','CeTD_SA3','DDOR_D','TPC_GKA']
        df_result = df_result.reindex(columns=column_name)
       
    
    elif dataset == "NON_AMP":
        df_ATC = namp_atc(df)
        df_DDOR = namp_DDOR(df)
        df_RRI = namp_RAAC(df)
        df_SER = namp_SE_residue_level(df)
        df_SEP= namp_sep_wp(df)
        df_CTD = namp_ctd(df)
        df_PAAC = namp_paac_1(df,3)
        df_APAAC = namp_apaac_1(df,3)
        df_QSO = namp_qos(df,3)
        df_SOC = namp_soc(df,3)
        df_result = pd.concat([df_ATC,df_DDOR,df_RRI,df_SER,df_SEP,df_CTD,df_PAAC,df_APAAC,df_QSO,df_SOC],axis = 1)
        column_names = ["ATC_H","APAAC3_HL_lam2","Negatively charged","Positively charged","ATC_N","CeTD_23_PO","PAAC3_lam2","CeTD_CH3","RRI_K",'Solvent Accesibilty(Intermediate)',"CeTD_CH1","DDOR_D","DDOR_E","QSO3_SC_K",'Acidicity','Basicity',
        "CeTD_12_VW",'Neutral',"SOC3_SC2","CeTD_31_VW","Polarity","CeTD_100_p_SS2","CeTD_50_p_PO1","CeTD_100_p_HB3","SER_K"]
        df_result = df_result.reindex(columns=column_names)
    return df_result
def mrmr_Predict_ALL_Features(dataset,method,sequence):
    
    
    if method == "predict":
        Fasta_df = [i[1].upper() for i in sequence ]
        Seq = [i[0] for i in sequence]
        df = pd.DataFrame(Fasta_df)
        df_result = mrmr_dataset_check(df,dataset)
        
        
    elif method == "protein_scan":
        Fasta_df = [i[1].upper() for i in sequence ]
        Seq = [i[0] for i in sequence]
        df = pd.DataFrame(Fasta_df)
        df_result = mrmr_dataset_check(df,dataset)
        
    
    elif method == "design": 
        Fasta_df = [i[1].upper() for i in sequence ]
        Seq = [i[0] for i in sequence]
        df = pd.DataFrame(Fasta_df)
        df_result = mrmr_dataset_check(df,dataset)
    return Fasta_df,Seq,df_result