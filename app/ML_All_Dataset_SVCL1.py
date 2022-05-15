import pandas as pd
import joblib, os
from .svc_all_features import Predict_ALL_Features
import numpy as np


def Model_pred(dataset,method,fastafile,model,threshold):
    if dataset == "Antiviral":
        if model == "SVM":
            model_file = joblib.load(os.path.join("SVC_model", "Antivral_Models", "AV_SVM_model.sav"))
        elif model == "DT":
            model_file = joblib.load(os.path.join("SVC_model", "Antivral_Models","AV_Decision_Tree_model.sav"))
        elif model == "RF":
            model_file = joblib.load(os.path.join("SVC_model", "Antivral_Models", "AV_Random_Forest_model.sav"))
        elif model == "Logistic_Regression":
            model_file = joblib.load(os.path.join("SVC_model", "Antivral_Models", "AV_Logistic_Regression_model.sav"))
        elif model == "XgBoost":
            model_file = joblib.load(os.path.join("SVC_model" , "Antivral_Models", "AV_XgBoost_model.sav"))
        
    elif dataset == "Antibacterial":
        if model == "SVM":
            model_file = joblib.load(os.path.join("SVC_model","Antibacterial_Models", "AB_SVM_model.sav"))
        elif model == "DT":
            model_file = joblib.load(os.path.join("SVC_model", "Antibacterial_Models", "AB_Decision_Tree_model.sav"))
        elif model == "RF":
            model_file = joblib.load(os.path.join("SVC_model", "Antibacterial_Models", "AB_Random_Forest_model.sav"))
        elif model == "Logistic_Regression":
            model_file = joblib.load(os.path.join("SVC_model","Antibacterial_Models","AB_Logistic_Regression_model.sav"))
        elif model == "XgBoost":
            model_file = joblib.load(os.path.join("SVC_model","Antibacterial_Models","AB_XgBoost_model.sav"))
        
    elif dataset == "Antifungal":
        if model == "SVM":
            model_file = joblib.load(os.path.join("SVC_model","Antifungal_Models", "AF_SVM_model.sav"))
        elif model == "DT":
            model_file = joblib.load(os.path.join("SVC_model","Antifungal_Models", "AF_Decision_Tree_model.sav"))
        elif model == "RF":
            model_file = joblib.load(os.path.join("SVC_model","Antifungal_Models", "AF_Random_Forest_model.sav"))
        elif model == "Logistic_Regression":
            model_file = joblib.load(os.path.join("SVC_model","Antifungal_Models", "AF_Logistic_Regression_model.sav"))
        elif model == "XgBoost":
            model_file = joblib.load(os.path.join("SVC_model","Antifungal_Models", "AF_XgBoost_model.sav"))
        
    elif dataset == "AMP":
        if model == "SVM":
            model_file = joblib.load(os.path.join("SVC_model","AMP_Models", "AMP_SVM_model.sav"))
        elif model == "DT":
            model_file = joblib.load(os.path.join("SVC_model","AMP_Models", "AMP_Decision_Tree_model.sav"))
        elif model == "RF":
            model_file = joblib.load(os.path.join("SVC_model","AMP_Models", "AMP_Random_Forest_model.sav"))
        elif model == "Logistic_Regression":
            model_file = joblib.load(os.path.join("SVC_model","AMP_Models", "AMP_Logistic_Regression_model.sav"))
        elif model == "XgBoost":
            model_file = joblib.load(os.path.join("SVC_model","AMP_Models", "AMP_XgBoost_model.sav"))
        
    elif dataset == "NON_AMP":
        if model == "SVM":
            model_file = joblib.load(os.path.join("SVC_model","d_NonAMP_Models", "d_Non_AMP_SVM_model.sav"))
        elif model == "DT":
            model_file = joblib.load(os.path.join("SVC_model","d_NonAMP_Models", "d_Non_AMP_Decision_Tree_model.sav"))
        elif model == "RF":
            model_file = joblib.load(os.path.join("SVC_model","d_NonAMP_Models", "d_Non_AMP_Random_Forest_model.sav"))
        elif model == "Logistic_Regression":
            model_file = joblib.load(os.path.join("SVC_model","d_NonAMP_Models", "d_Non_AMP_Logistic_Regression_model.sav"))
        elif model == "XgBoost":
            model_file = joblib.load(os.path.join("SVC_model","d_NonAMP_Models", "d_Non_AMP_XgBoost_model.sav"))

    Fasta_Seq,ID,test_full_features = Predict_ALL_Features(dataset, method, fastafile)
    
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
  