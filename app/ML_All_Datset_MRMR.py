import pandas as pd
import joblib, os
from .mrmr_all_features import mrmr_Predict_ALL_Features
import numpy as np

def MRMR_Model_pred(dataset,method,fastafile,model,threshold):
    if dataset == "Antiviral":
        if model == "SVM":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiViral", "MRMR_AV_SVM_model.sav"))
        elif model == "DT":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiViral", "MRMR_AV_Decision_Tree_model.sav"))
        elif model == "RF":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiViral", "MRMR_AV_Random_Forest_model.sav"))
        elif model == "Logistic_Regression":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiViral", "MRMR_AV_Logistic_Regression_model.sav"))
        elif model == "XgBoost":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiViral", "MRMR_AV_XgBoost_model.sav"))
        
    elif dataset == "Antibacterial":
        if model == "SVM":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiBacterial", "MRMR_AB_SVM_model.sav"))
        elif model == "DT":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiBacterial", "MRMR_AB_Decision_Tree_model.sav"))
        elif model == "RF":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiBacterial", "MRMR_AB_Random_Forest_model.sav"))
        elif model == "Logistic_Regression":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiBacterial", "MRMR_AB_Logistic_Regression_model.sav"))
        elif model == "XgBoost":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiBacterial", "MRMR_AB_XgBoost_model.sav"))
        
    elif dataset == "Antifungal":
        if model == "SVM":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiFungal", "MRMR_AF_SVM_model.sav"))
        elif model == "DT":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiFungal", "MRMR_AF_Decision_Tree_model.sav"))
        elif model == "RF":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiFungal", "MRMR_AF_Random_Forest_model.sav"))
        elif model == "Logistic_Regression":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiFungal", "MRMR_AF_Logistic_Regression_model.sav"))
        elif model == "XgBoost":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AntiFungal", "MRMR_AF_XgBoost_model.sav"))
        
    elif dataset == "AMP":
        if model == "SVM":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AMP", "MRMR_AMP_SVM_model.sav"))
        elif model == "DT":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AMP", "MRMR_AMP_Decision_Tree_model.sav"))
        elif model == "RF":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AMP", "MRMR_AMP_Random_Forest_model.sav"))
        elif model == "Logistic_Regression":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AMP", "MRMR_AMP_Logistic_Regression_model.sav"))
        elif model == "XgBoost":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_AMP", "MRMR_AMP_XgBoost_model.sav"))
        
    elif dataset == "NON_AMP":
        if model == "SVM":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_NonAMP", "MRMR_NonAMP_SVM_model.sav"))
        elif model == "DT":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_NonAMP", "MRMR_NonAMP_Decision_Tree_model.sav"))
        elif model == "RF":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_NonAMP", "MRMR_NonAMP_Random_Forest_model.sav"))
        elif model == "Logistic_Regression":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_NonAMP", "MRMR_NonAMP_Logistic_Regression_model.sav"))
        elif model == "XgBoost":
            model_file = joblib.load(os.path.join("MRMR_model", "MRMR_NonAMP", "MRMR_NonAMP_XgBoost_model.sav"))

    Fasta_Seq,ID,test_full_features = mrmr_Predict_ALL_Features(dataset, method, fastafile)
    
    Result = []
    
    for index,row in test_full_features.iterrows():
        row = pd.DataFrame(row)
        row_T = row.transpose()
        row_T = np.array(row_T)
        probability = model_file.predict_proba(row_T)
        probability = np.round(probability,4)
        if probability[0][1] >= threshold:
            Result.append([ID[index],Fasta_Seq[index],"Antiprotozoal Activity",probability[0][1]])
        else:
            Result.append([ID[index],Fasta_Seq[index],"No Antiprotozoal Activity",probability[0][1]])
    return Result
