import pickle
import pybel as pb
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors


class PredictBpSol(object):
    '''
    Provide a user interface for predicting boiling points of alkanes.
    Input -> smiles string, Output -> boiling point
    '''
    def __init__(self, smiles):
        self.smiles = smiles

    def predict(self, path_to_model):
        descr_names = [i[0] for i in Descriptors._descList if i[0] != 'ExactMolWt']
        m = Chem.MolFromSmiles(str(self.smiles))
        calc = MoleculeDescriptors.MolecularDescriptorCalculator(descr_names)
        descr = [i for i in calc.CalcDescriptors(m)]
        with open(path_to_model, 'rb') as f:
            self.model = pickle.load(f)
        self.predicted = self.model.predict(descr)
        return '{:.2f}'.format(self.predicted[0])


all_toxi = ['a[N;X2]=O', 'CO[N;X2]=O', 'N[N;X2]=O', 'O1[c,C]-[c,C]1',
            'C1NC1', 'N=[N+]=[N-]', 'C=[N+]=[N-]', 'N=N-N',
            'cN!@;=[N;X3]', '[OH]O', '[OH][N;X2]=C', '[c,C]OO[c,C]',
            'C[NH][NH]C', '[OH]Na', '[NH2]Na', '[OH][N;X2]=[N;X2]',
            '[Cl,Br,I]C', '[Cl,Br,I]C=O', '[N,S]!@[C;X4]!@[CH2][Cl,Br,I]',
            'SC[Cl]', '[Cl,Br,I]!@[C;X4]!@[C;X4]O', '[Cl]C([X1])=C[X1]',
            '[Cl,Br,I]C(([F,Cl,Br,I])[X1])C=C',
            '[cH]1[cH]ccc2c1c3c(cc2)cc[cH][cH]3',
            '[Cl,Br,I]C(([F,Cl,Br,I])[X1])C(=O)[c,C]',
            '[cH]1ccc2c1[cH][cH]c3c2ccc[cH]3', '[C,c]OS((=O)=O)O!@[c,C]',
            '[c,C]S((=O)=O!@[c,C])', 'O=N(~O)N', '[N;v4]#N', 'O=C1CCO1',
            '[CH]=[CH]O', 'aN([OH])', 'aN(O*=O)', 'aN(C=O)[CH3]',
            'aN[CH3]', 'aN([#1])', '[NH;!R][NH;!R]a', '[CH3][NH]a',
            'a13~a~a~a~a2~a1~a(~a~a~a~3)~a~a~a~2',
            'a1~a~a~a2~a~1~a~a3~a(~a~2)~a~a~a~3',
            'a1~a~a~a2~a~1~a~a~a3~a~2~a~a~a~3',
            'a1~a~a~a~2a~a~1~a3~a(~a~2)~a~a~a~a~3',
            'a1~a~a~a~2a~a~1~a~a3~a(~a~2)~a~a~a~3',
            'a1~a~a~a~2a~a~1~a~a3~a(~a~2)~a~a~a~a~3',
            'a1~a~a~a~2a~a~1~a~a~a3~a~2~a~a~a~3',
            'a1~a~a~a~2a~a~1~a~a~a3~a~2~a~a~a~a~3',
            'a13~a~a~a~2a~a1~a(~a~a~a~3)~a~a~2']

d = {'O=N(~O)a': ['O=N(O)c(aS(=O)=O)', 'O=N(O)c(aaS(=O)=O)',
                  'O=N(O)c(aaaS(=O)=O)', 'O=N(O)c(aC((F)F)F)',
                  'O=N(O)c(aaC((F)F)F)', 'O=N(O)c(aaaC((F)F)F)'],
     'a[NH2]': ['[NH2]a(a(C((F)F)F))', '[NH2]a(a(S(=O)=O))',
                '[NH2]a(a(C(=O)O))', '[NH2]a(aa(C((F)F)F))',
                '[NH2]a(aa(S(=O)=O))' '[NH2]a(aa(C(=O)O))',
                '[NH2]a(aaa(C((F)F)F))', '[NH2]a(aaa(S(=O)=O))',
                '[NH2]a(aaa(C(=O)O))'],
     'c[N;X2]!@;=[N;X2]c': ['[N;X2](acS((=O)=O))=[N;X2](acS((=O)=O))',
                            '[N;X2](aacS((=O)=O))=[N;X2](aacS((=O)=O))',
                            '[N;X2](aaacS((=O)=O))=[N;X2](aaacS((=O)=O))',
                            '[N;X2](aaaacS((=O)=O))=[N;X2](aaaacS((=O)=O))'],
     '[OH,NN2][N,O]': ['O=N(O)[O-]'],
     '[OH]N': ['[OH]Na', '[OH][N;X2]=*'],
     '[Cl,Br,I][C;X4]': ['[Cl,Br,I][C;X4][F,Cl,Br,I]',
                         '[Cl,Br,I]C((C)C)C'],
     '[Cl,Br,I][CH][CH3]': ['[Cl,Br,I][C]([Cl,Br,I,F])[CH3]'],
     'O=[CH]C=C': ['O=[CH]C([N,O,S])=C', 'O=[CH]C=C[N,O,S]',
                   'O=[CH]C=Ca'],
     'O=[CH]C=O': ['O=[CH]C([N,O,S])=C', 'O=[CH]C=C[N,O,S]',
                   'O=[CH]C=Ca'],
     '[CH3][NH]a': ['[CH3][NH]a(a(C((F)F)F))', '[CH3][NH]a(a(S=O))',
                    '[CH3][NH]a(a(C(=O)O))', '[CH3][NH]a(aa(C((F)F)F))',
                    '[CH3][NH]a(aa(S=O))', '[CH3][NH]a(aa(C(=O)O))',
                    '[CH3][NH]a(aaa(C((F)F)F))', '[CH3][NH]a(aaa(S=O))',
                    '[CH3][NH]a(aaa(C(=O)O))']
     }


class Mutagenicity(object):
    '''Mutagenicity classifier'''
    def __init__(self, smiles):
        self.smiles = smiles

    def predict(self, path_to_model):
        descr_names = [i[0] for i in Descriptors._descList if i[0] != 'ExactMolWt']
        m = Chem.MolFromSmiles(str(self.smiles))
        calc = MoleculeDescriptors.MolecularDescriptorCalculator(descr_names)
        descr = [i for i in calc.CalcDescriptors(m)]
        with open(path_to_model, 'rb') as f:
            self.model = pickle.load(f)
        self.predicted = self.model.predict(descr)
        return self.predicted[0]

    def unique_toxi(self):
        mol_list = []
        for smi in [self.smiles]:
            l = []
            mol = pb.readstring("smi", str(smi))
            for toxi in all_toxi:
                smarts = pb.Smarts(toxi)
                if len(smarts.findall(mol)) > 0:
                    l.append(1)
                else:
                    l.append(0)
            mol_list.append(sum(l))
        return mol_list

    def dict_toxi(self):
        mol_list = []
        for smi in [self.smiles]:
            l = []
            mol = pb.readstring("smi", str(smi))
            for k, v in d.iteritems():
                k_smarts = pb.Smarts(k)
                n = len(k_smarts.findall(mol))
                if n == 0:
                    l.append(0)
                else:
                    for each in v:
                        d_list = []
                        v_smarts = pb.Smarts(each)
                        d_list.append(len(v_smarts.findall(mol)))
                        if n > sum(d_list):
                            l.append(1)
                        elif n == sum(d_list):
                            l.append(0)
            mol_list.append(sum(l))
        return mol_list

    def merge_toxi(self):
        all_toxi_labels = self.unique_toxi()[0] + self.dict_toxi()[0]
        if all_toxi_labels == 0:
            return 0
        else:
            return 1

    def vote_for_spec(self, path_to_model):
        if self.predict(path_to_model) == 1:
            return 1
        else:
            return self.merge_toxi()


def get_image(smiles, path_to_image, legend):
    smi = [str(smiles)]
    mol = [AllChem.MolFromSmiles(sm) for sm in smi]
    img = Draw.MolsToGridImage(mol, molsPerRow=1,
                               subImgSize=(200, 200),
                               legends=[legend])
    img.save(path_to_image)
