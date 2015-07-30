from flask import Flask, render_template, request, jsonify
from forms import forms
from chemistry.chemfig import smiles_mol_to_chemfig, get_name, update_chemfig
import re
from models import ml_models
import pybel
import os


app = Flask(__name__)
app.config.from_object('config.DevelopmentConfig')

#################
# Mol_2_chemfig #
#################

MOL_FILE = 'static/files/molecule.mol'


@app.route('/mol_2_chemfig')
def home():
    form = forms.HomePageForm()
    global text_area
    text_area = form.smiles_mol.data
    checkboxes = forms.get_checkboxes_data()
    menu_links = forms.get_menu_links()
    return render_template('home.html', menu_links=menu_links, form=form,
                           checkboxes=checkboxes, chkbx_class="check_box",
                           label_class="checkboxtext",
                           pdflink='/static/files/welcome.pdf')


@app.route('/mol_2_chemfig/links')
def links():
    return render_template("links.html")


@app.route('/mol_2_chemfig/about')
def about():
    return render_template("about.html")


# Keeps track of the last content in the text area
last_data = ["text"]


@app.route("/mol_2_chemfig/_get_smiles")
def get_smile():
    chemical = request.args.get('chemical')
    name = get_name(chemical)
    return jsonify(smiles=name)


def get_options():
    smiles_mol = request.args.get("smiles_mol")
    last_data[0] = smiles_mol
    check = request.args.getlist('check')
    check = ' '.join(check[0].split(','))
    angle = request.args.get('angle')
    angle = " -a " + angle
    lst = check + angle
    hydrogens = request.args.get("hydrogens")
    return smiles_mol, lst, hydrogens


def display(smiles_mol, lst, hydrogens, update=False):
    if len(smiles_mol) < 200:
        chemfig, pdflink = smiles_mol_to_chemfig(lst, '-i direct',
                                                 "-y {}".format(hydrogens),
                                                 smiles_mol)
    else:
        if not update:
            with open(MOL_FILE, 'w') as f:
                f.write(smiles_mol)
        chemfig, pdflink = smiles_mol_to_chemfig(lst,
                                                 "-y {}".format(hydrogens),
                                                 MOL_FILE)
    return chemfig, pdflink


@app.route("/mol_2_chemfig/smiles_to_chemfig")
def smiles_to_chemfig():
    smiles_mol, lst, hydrogens = get_options()
    chemfig, pdflink = display(smiles_mol, lst, hydrogens)
    s = re.sub('\%.+', '', chemfig)

    def inline():
        my_string = ''
        for i in s.split():
            my_string = my_string + i
        return my_string

    if '-r' in lst:
        chemfig = s
    elif '-in' in lst:
        chemfig = inline()
    return jsonify(chem_fig=chemfig, pdf_link=pdflink)


@app.route("/mol_2_chemfig/update")
def check_update():
    smiles_mol, lst, hydrogens = get_options()
    chemfig, pdflink = display(smiles_mol, lst, hydrogens, update=True)
    if '-r' in lst:
        chemfig = re.sub('\%.+', '', chemfig)
    return jsonify(chem_fig=chemfig, pdf_link=pdflink)


@app.route("/mol_2_chemfig/update_chemfig")
def chemfig_update():
    smiles_mol = request.args.get("smiles_mol")
    pdflink = update_chemfig(smiles_mol)
    return jsonify(pdf_link=pdflink)


####################
# Machine learning #
####################

MOL_FILE = 'static/files/molecule2.mol'
BP = 'static/models/bp_model.pickle'
SOL = 'static/models/sol_model.pickle'
MP = 'static/models/mp_model.pickle'
MUT = 'static/models/mutagenicity_model.pickle'
PATH_TO_IMAGE = 'static/img/ml_app'


@app.route('/ml_app')
def home_ml():
    Form = forms.HomePageFormML()
    return render_template('home_ml.html',
                           form=Form,
                           image=' ',
                           units='',
                           symbol='')


@app.route('/ml_app/about')
def about_ml():
    return render_template('about_ml.html')


@app.route('/ml_app/get_smiles')
def get_smiles_ml():
    smiles_mol = request.args.get('smiles_mol')
    with open(MOL_FILE, 'w') as f:
        f.write(smiles_mol)
    smiles = pybel.readfile('mol', MOL_FILE).next()
    smi = smiles.write('smi').split('\t')[0]
    # can_smi = pybel.readstring('smi', smi).write('can')
    return jsonify(smiles=smi)


def del_image_ml():
    files = os.listdir(PATH_TO_IMAGE)
    img = [i for i in files if i.endswith('png')]
    if img:
        for i in img:
            os.remove(PATH_TO_IMAGE + '/' + i)


@app.route('/ml_app/predicted')
def predict_value():
    del_image_ml()
    img_file = PATH_TO_IMAGE + '/{}.png'
    smiles = request.args.get('smiles')
    ml_options = request.args.get('ml_options')
    model = ml_models.PredictBpSol(smiles)
    if ml_options == 'bp':
        value = float(model.predict(BP))
        if -200.0 < value < 450.0:
            result = 'Predicted value: {}'.format(value)
            units = 'Units: Degrees Celsius'
        else:
            units = "Value can't be generated"
            result = 'Please use alkanes'
    elif ml_options == 'sol':
        result = 'Predicted value: {}'.format(model.predict(SOL))
        units = 'Units: log(mol/L)'
    elif ml_options == 'toxi':
        model_toxi = ml_models.Mutagenicity(smiles)
        units = ''
        if model_toxi.vote_for_spec(MUT) == 0:
            result = 'Nonmutagenic'
        else:
            result = 'Mutagenic'
    elif ml_options == 'mp':
        result = 'Predicted value: {}'.format(model.predict(MP))
        units = 'Units: Degrees Celsius'
    else:
        result = 'Sorry'
    img_name = smiles.replace('/', '').replace('\\', '').replace('#', '')
    ml_models.get_image(smiles, img_file.format(result+img_name), result)
    return jsonify(image=img_file.format(result+img_name),
                   units=units)

if __name__ == '__main__':
    app.run()
