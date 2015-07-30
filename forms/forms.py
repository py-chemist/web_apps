from flask.ext.wtf import Form
from wtforms import StringField, TextAreaField
from wtforms import SubmitField, BooleanField, SelectField

#################
# Mol_2_chemfig #
#################


def get_checkboxes_data():
    '''Return title, values and text of checkboxes'''

    title = ["Wrap the code into \chemfig{...}",
             "Remove %",
             "Dispaly chemfig code inline",
             "Assign number to each atom except hydrogen",
             "Show nicer double and triple bonds",
             "Show circle in aromatic compounds instead of double bonds",
             "Show carbon atoms as elements",
             "Show methyl group as elements",
             "Flip the structure horizontally",
             "Flop the strcture vertically"]
    value = ["-w", "-r", "-in", "-n", "-f", "-o", "-c", "-m", "-p", "-q"]
    text = ["chemfig", "remove %", "inline-mode",
            "atom-numbers", "fancy bonds", "aromatic",
            "show carbon", "show methyl", "flip", "flop"]
    return zip(title, value, text)


def get_menu_links():
    return ["Home", "Links", "About"]


def hydrogens_options():
    return [("keep", "keep"), ("add", "add"), ("delete", "delete")]


class HomePageForm(Form):
    db_search = StringField()
    submit = SubmitField('Search')
    smiles_mol = TextAreaField()
    action = SubmitField("Convert")
    check = BooleanField()
    angle = StringField()
    hydrogens = SelectField('hydrohens', choices=hydrogens_options())
    update = SubmitField('Apply')
    reset = SubmitField('Reset')


####################
# Machine learning #
####################


def choose_ml():
    return [('sol', 'Aqueous solubility of organic compounds'),
            ('toxi', 'Mutagenicity of organic compounds'),
            ('mp', 'Melting point of organic compounds '),
            ('bp', 'Boiling point of alkanes')
            ]


class HomePageFormML(Form):
    ml_options = SelectField('ml_options', choices=choose_ml())
    smiles = StringField()
    run = SubmitField('Run')
