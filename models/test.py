import ml_models

model = ml_models.Mutagenicity('NC(CCC(=O)NC(CSC(=O)NCCCl)C(=O)NCC(=O)O)C(=O)O')

#model.predict('/home/py-chemist/Projects/Machine_Learning/Mutagenes/mutagen_model.pickle')

#print model.unique_toxi()

#print model.dict_toxi()

print model.vote_for_spec('/home/py-chemist/Projects/Machine_Learning/Mutagenes/mutagen_model.pickle')
