def result_message(mess):
    if mess == "predict":
        return "This page display the results of the Predict Module for the query sequences submitted by the user. The table below shows the details of the query peptide sequences, with the first column displaying the serial number of the input peptide sequence, the second column displaying the ID of the submitted input sequence, third column representing the sequence of the submitted peptide, fourth column displaying the probability score of the given peptide being antiprotozoal as obtained by the user defined negative dataset, feature selection and machine learning model. The last column display the prediction whether the peptide might have antiprotozoal activity or not determined by the condition whether the probability score is greater or less than the user defined threshold."
    elif mess == "protein":
        return "This page display the results of the protein scan module of the query protein. The table below shows the details of the query protein sequences, with the first column displaying the serial number, the second column displaying the ID of the peptide fragment, the third column displaying the sequence of the generated peptide, fourth column displaying the antiprotozoal probability as provided by the machine learning algorithm based on various parameters (negative dataset, feature selection and prediction model) selected by the user. The last column display the prediction whether the peptide might have antiprotozoal activity or not determined by the condition whether the probability score is greater or less than the user defined threshold."

    elif mess == "design":
        return "This page display the results of the design module in the single substitution mutants of the query peptide. The table below shows the details of the query peptide sequences, with the first column displaying the serial number, the second column displaying the ID of the mutant sequence, the third column displaying the sequence of the mutant generated, fourth column displaying the probability provided by the machine learning algorithm based on the negative dataset, feature selection and prediction model chosen by the user. Last column display the prediction whether the peptide and its mutant might have antiprotozoal activity or not determined by the condition whether the probability score is greater or less than the user defined threshold."


def result_title(mess):
    if mess == "predict":
        return "Result of Predict Module"
    elif mess == "protein":
        return "Result of Protein Scan Module"

    elif mess == "design":
        return "Result of Design Module"

def error_message(mess):
    if mess == "404":
        return "oops!! Page not found 404!!"
    elif mess == "500":
        return "oops!! Something Went Wrong!!"

def error_title(mess):
    if mess =="404":
        return "404"
    elif mess =="500":
        return "500"