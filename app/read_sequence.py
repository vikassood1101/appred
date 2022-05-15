import os 
from flask import flash
import re
import string


def predict_validation(sequence):
    sequence = sequence.strip()
    sequence = sequence.split("\r\n")
    invalid_char = set(string.punctuation)

    for i in range(0,len(sequence),2):
        if not re.search('>',sequence[i]):
            return flash("The uploaded file/input sequence(s) is not in fasta format","danger")
    
    records = [i.rstrip("\n") for i in sequence]
    for fasta in range(1,len(records),2):
        if any(char in invalid_char for char in records[fasta]):
            return flash("Your sequence(s) contains special characters","danger")
            
    for fasta in range(1,len(records),2):
        if " " in records[fasta]:
            return flash("Your sequence(s) contains white space(s)","danger")

    for fasta in range(1,len(records),2):
       if records[fasta].isalpha() and  re.search('[^acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY]',records[fasta]):
            return flash("Your sequence(s) contains non-natural amino acid(s)","danger")

    for fasta in range(1,len(records),2):
        if not records[fasta].isalpha() :
            return flash("Your sequence(s) contains digits. Please enter only alphabetical characters","danger")         

    for fasta in range(1,len(records),2):
        if len(records[fasta]) < 8 or len(records[fasta]) > 100:
            return flash("The length of the peptide sequence(s) is not in the range of 8-100 amino acids","danger")
    records_tuple = [(records[x], records[x+1]) for x in range(0,len(records),2)]
    
    return records_tuple

##################Protein_Scan###############

def protein_validation(sequence,k):
    sequence = sequence.strip()
    sequence = sequence.split("\r\n")
    if len(sequence)>1:
        return flash("Please provide only one protein sequence in single line","danger")
    sequence = sequence[0]
    invalid_char = set(string.punctuation)
    if any(char in invalid_char for char in sequence):
        return flash('Your sequence contains special character(s)',"danger")
        
    if " " in sequence:
        return flash("Your sequence contains white space(s)","danger")
   
    if len(sequence) < 8 or len(sequence) > 100:
        return flash("The length of the peptide sequence is not in the range of 8-100 amino acids","danger")   
    
    if sequence.isalpha() and  re.search('[^acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY]',sequence):
        return flash("Your sequence contains non-natural amino acid(s)","danger")

    if not sequence.isalpha():
        return flash("Your Sequence contains digits. Please enter only alphabetical characters","danger")     


    if len(sequence) < k:
        return flash("The length of peptide sequence is less than the peptide fragment length","danger")

     
    
    res = [sequence[idx:idx + k] for idx in range(len(sequence) - k + 1)]
    Seq = [f"Seq {x+1}" for x in range(len(res))]
    return list(zip(Seq, res)) 

##################### DESIGN ##################

def motif_validation(sequence):
    sequence = sequence.strip()
    list_Mutant = []
    characters = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    invalid_char = set(string.punctuation)
    sequence = sequence.split("\r\n")
    if len(sequence)>1:
        return flash("Please provide only one protein sequence in single line, as shown in the example.","danger")
    sequence = sequence[0]
    if any(char in invalid_char for char in sequence):
        return flash('Your sequence contains special character(s)',"danger")

    if " " in sequence:
        return flash("Your sequence contain white space(s)","danger")
   
    if len(sequence) < 8 or len(sequence) > 100:
        return flash("The length of the peptide sequence is not in the range of 8-100 amino acids","danger")

    if sequence.isalpha() and  re.search('[^acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY]',sequence):
        return flash("Your sequence contains non-natural amino acid(s)","danger")

    if not sequence.isalpha():
        return flash("Your sequence contains digits. Please enter only alphabetical characters","danger")    

    sequence = sequence.upper()    
    for i in range(len(sequence)):
        for j in characters:
            if i != j:
                #sample_str = sequence[0: i] + j + sequence[i + 1:]
                list_Mutant.append(sequence[0: i] + j + sequence[i + 1:])
            else:
                list_Mutant.append(sequence)
    Seq = [f"Seq {x+1}" for x in range(len(list_Mutant))]
    return list(zip(Seq, list_Mutant)) 
