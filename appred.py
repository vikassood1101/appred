
from flask import Flask ,render_template, redirect, request,flash
import logging
from app.read_sequence import predict_validation, protein_validation,motif_validation
from app.ML_All_Dataset_SVCL1 import *
from app.svc_all_features import *
from app.ML_All_Datset_MRMR import *
from app.mrmr_all_features import *
from app.message import result_message,result_title,error_message,error_title

app = Flask(__name__)

app.config.update({'SECRET_KEY' : 'FQ6NepofYh2aBgK94nbpwHkyBEs2b5NF3TF3z3XlbMxLezVQHe'})

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename= os.path.join(app.static_folder, 'logger', 'appred.log'),
                    filemode='a')

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/predict", methods = ["GET","POST"])
def predict():
    import sys
    logging.info(sys.path[1])
    logging.info(os.getcwd())
    logging.info(os.curdir)
    ROOT_DIR = os.path.abspath(os.curdir)
    logging.info(ROOT_DIR)
    if request.method =="POST":
        form = request.form
        sequence = form["sequence"]
        negative_data = form["negdata"]
        feature_selection = form["featureSelection"]
        modal_selection  = form["ModalSelection"]
        threshold = float(form["threshold"])
        if len(sequence) == 0:
            flash("sequence can't be empty","danger")
            return render_template("predict.html")
        file = predict_validation(sequence)
        if file == None:
            return render_template("predict.html")
        if feature_selection =="SVC":
            result = Model_pred(negative_data,"predict",file,modal_selection,threshold)
            return render_template("results.html",tables = result,result_message=result_message("predict"),result_title=result_title("predict"))
        elif feature_selection =="MRMR":
            result = MRMR_Model_pred(negative_data,"predict",file,modal_selection,threshold)
            return render_template("results.html",tables = result,result_message=result_message("predict"),result_title=result_title("predict"))
    return render_template("predict.html")

@app.route("/proteinScan",methods=["GET","POST"])
def protien():
    if request.method == "POST":
        form = request.form
        sequence = form["sequence"]
        negative_data = form["negdata"]
        feature_selection = form["featureSelection"]
        modal_selection  = form["ModalSelection"]
        threshold = float(form["threshold"])
        k = int(form["slider_value"])
        if len(sequence) == 0:
            flash("sequence can't be empty","danger")
            return render_template("protein.html")
        file = protein_validation(sequence, k)
        if file == None:
            return render_template("protein.html")
        if feature_selection =="SVC":
            result = Model_pred(negative_data,"predict",file,modal_selection,threshold)
            return render_template("results.html",tables = result,result_message=result_message("protein"),result_title=result_title("protein"))
        elif feature_selection =="MRMR":
            result = MRMR_Model_pred(negative_data,"predict",file,modal_selection,threshold)
            return render_template("results.html",tables = result,result_message=result_message("protein"),result_title=result_title("protein"))
    return render_template("protein.html")

@app.route("/dataset")
def dataset():
    return render_template("dataset.html")

@app.route("/help")
def help():
    return render_template("help.html")

@app.route("/contact")
def contact():
    return render_template("contact.html")

@app.route("/design",methods=["GET","POST"])
def design():
    if request.method =="POST":
        form = request.form
        sequence = form["sequence"]
        negative_data = form["negdata"]
        feature_selection = form["featureSelection"]
        modal_selection  = form["ModalSelection"]
        threshold = float(form["threshold"])
        if len(sequence) == 0:
            flash("sequence can't be empty","danger")
            return render_template("design.html")
        file = motif_validation(sequence)
        if file == None:
            return render_template("design.html")
        if feature_selection =="SVC":
            result = Model_pred(negative_data,"predict",file,modal_selection,threshold)
            return render_template("results.html",tables = result,result_message=result_message("design"),result_title=result_title("design"))
        elif feature_selection =="MRMR":
            result = MRMR_Model_pred(negative_data,"predict",file,modal_selection,threshold)
            return render_template("results.html",tables = result,result_message=result_message("design"),result_title=result_title("design"))
    return render_template("design.html")

@app.errorhandler(404)
def page_not_found(e):
    return render_template("error.html",error_title=error_title("404"),error_message=error_message("404"))

@app.errorhandler(500)
def page_not_found(e):
    return render_template("error.html",error_title=error_title("500"),error_message=error_message("500"))

if __name__== "__main__":
    app.run(debug=True)