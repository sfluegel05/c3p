"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
for carbon in connected_carbons:
            if not carbon.HasSubstructMatch(Chem.MolFromSmarts('[CX4][OX2H]')):
              current_atom = carbon
              fatty_acid_carbon_chain.append(current_atom.GetIdx())
              break