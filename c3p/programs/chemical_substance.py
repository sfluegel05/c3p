"""
Classifies: CHEBI:59999 chemical substance
"""
from rdkit import Chem

def is_chemical_substance(smiles: str):
    """
    Determines if a molecule is a chemical substance (portion of matter of constant composition,
    composed of molecular entities of the same type or of different types).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chemical substance, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule is composed of molecular entities
    if mol.GetNumAtoms() == 0:
        return False, "No atoms found in the molecule"
    
    # For simplicity, assume any valid molecule with atoms is a chemical substance
    return True, "The molecule is a chemical substance"

# Example usage
print(is_chemical_substance("CCCCCCCCCCCCCCCCC/C=C\CCCCCCCCCCCCC/C=C/CCCCCCCCCCCCCCCCCC(C(C(=O)O)CCCCCCCCCCCCCCCCCCCCCC)O"))  # alpha-Semegma mycolic acid
print(is_chemical_substance("CCCCCC\C=C/C\C=C/CCCCCCCC(=O)OCC(O)CO"))  # rac-1-monolinoleoylglycerol
print(is_chemical_substance("CC(C)=CCC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\Cc1ccccc1O"))  # 2-all-trans-decaprenylphenol
print(is_chemical_substance("O.Cc1ccc(cc1)S(O)(=O)=O.NC1CCN(C1)c1nc2n(cc(C(O)=O)c(=O)c2cc1F)-c1ccc(F)cc1F"))  # tosufloxacin tosylate hydrate
print(is_chemical_substance("C(=S)([S-])N(CCCC)CCCC.C(=S)([S-])N(CCCC)CCCC.[Zn+2].C(=S)([S-])N(CC)CC.C(=S)([S-])N(CC)CC.[Zn+2].C(=N)(NC=1C=CC=CC1)NC2=CC=CC=C2"))  # carba mix
print(is_chemical_substance("[Ba++].[O-]C([O-])=O"))  # barium carbonate
print(is_chemical_substance("COc1ccc(cc1OC)C(=CC(=O)N1CCOCC1)c1ccc(Cl)cc1"))  # dimethomorph
print(is_chemical_substance("O(CCCCCCCCCCCCCCCCCCCCCCCCCCCC)C(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCCC"))  # Octacosyl triacontanoate
print(is_chemical_substance("CC(C)=CCC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\Cc1ccccc1O"))  # 2-all-trans-nonaprenylphenol
print(is_chemical_substance("C[C@@H](NC(=O)C(C#N)C(C)(C)C)c1ccc(Cl)cc1Cl"))  # diclocymet
print(is_chemical_substance("O(CCCCCCCC/C=C\C/C=C\C/C=C\CC)C(=O)CCCCCCC/C=C\C/C=C\CCCCC"))  # Linolenyl linoleate
print(is_chemical_substance("CCCCCCCC\C=C/CCCCCCCC(=O)OC(CO)COC(=O)CCCCCCC\C=C/C\C=C/CCCCC"))  # rac-1-linoleoyl-2-oleoylglycerol
print(is_chemical_substance("CCCCCCCCCCCCCCCCCC(=O)OCCCCCCCC"))  # octyl palmitate
print(is_chemical_substance("O(CCCCCCCCCCCCCCCCCCCCCC)C(=O)CCCCCCC/C=C\CCCCCCCC"))  # Behenyl oleate
print(is_chemical_substance("CCCCCCCCCCCCCCCCCCCOC(=O)CCCCCCCCCCCCCCC"))  # stearyl palmitate
print(is_chemical_substance("[Al](O[Si](O[Si](O[Al]=O)=O)=O)=O.O.O"))  # kaolin
print(is_chemical_substance("C1(=CC=C(C=C1)/C=C/C2=CC=C(C=C2)N)N.Cl.Cl"))  # 4,4'-Diaminostilbene dihydrochloride
print(is_chemical_substance("C1=CC=CC=C1C[C@@H](C(NCCCCCC(O)=O)=O)NC(=O)C2=CC(=C3C(=C2O)C(OC(C3)C)=O)Cl"))  # hapten OTAe
print(is_chemical_substance("O(CCCCCCCC/C=C\CCCC)C(=O)CCCCCCCCCCCCC"))  # Myristoleyl myristate
print(is_chemical_substance("C(=O)([C@@H](N)CCC(=O)N[C@@H](CC1(C(=C)C1)[H])C(O)=O)O"))  # hypoglycin B
print(is_chemical_substance("[Cs+].[O-][N+]([O-])=O"))  # caesium nitrate
print(is_chemical_substance("CCCCCCCC\C=C/CCCCCCCC(=O)OC(CO)COC(=O)CCCCCCC\C=C/C\C=C/C\C=C/CC"))  # rac-1-alpha-linolenoyl-2-oleoylglycerol
print(is_chemical_substance("C\C=C\C(=O)Oc1c([*])cc(cc1[N+]([O-])=O)[N+]([O-])=O"))  # dinocap-6
print(is_chemical_substance("O(CCCCCCCCCCCCCC)C(=O)CCCCCCC/C=C\C/C=C\CCCCC"))  # Myristyl linoleate
print(is_chemical_substance("CC(C)=CCC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\CC\C(C)=C\Cc1ccccc1O"))  # 2-all-trans-heptaprenylphenol
print(is_chemical_substance("*:c1c(:*)c2c(:*)c(:*)c3c(:*)c4c(:*)c(:*)c(:*)c5c(:*)c(:*)c6c(:*)c(c1:*)c2c3c6c45.*:c1c(:*)c2c(:*)c(:*)c3c(:*)c4c(:*)c(:*)c(:*)c5c(:*)c(:*)c6c(:*)c(c1:*)c2c3c6c45.*:c1c(:*)c2c(:*)c(:*)c3c(:*)c4c(:*)c(:*)c(:*)c5c(:*)c(:*)c6c(:*)c(c1:*)c2c3c6c45.*:c1c(:*)c2c(:*)c(:*)c3c(:*)c4c(:*)c(:*)c(:*)c5c(:*)c(:*)c6c(:*)c(c1:*)c2c3c6c45"))  # graphite
print(is_chemical_substance("O(CCCCCCCC/C=C\CCCC)C(=O)CCCCCCCCCCC"))  # Myristoleyl laurate
print(is_chemical_substance("O(CCCCCCCC/C=C\CCCCCCCC)C(=O)CCCCCCC/C=C\CCCCCC"))  # Oleyl palmitoleate
print(is_chemical_substance("CCN1OCC(NC(=O)C2=C(C)C=C(C=C2)C2=NOC(C2)(C2=CC(Cl)=C(F)C(Cl)=C2)C(F)(F)F)C1=O"))  # isocycloseram


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59999',
                          'name': 'chemical substance',
                          'definition': 'A chemical substance is a portion of '
                                        'matter of constant composition, '
                                        'composed of molecular entities of the '
                                        'same type or of different types.',
                          'parents': ['CHEBI:24431']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n"
              "(True, 'The molecule is a chemical substance')\n",
    'num_true_positives': 29,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.5918367346938775,
    'recall': 1.0,
    'f1': 0.7435897435897436,
    'accuracy': None}