"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside contains a pyrimidine base linked to a 2-deoxyribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pyrimidine base (generic) - a six-membered ring with 2 nitrogens and 4 carbons.
    # One nitrogen should be at position 1. The other nitrogen can be at position 3.
    # Accept any substituents on the ring.
    # Also include explicit H counts to accept protonated rings.

    pyrimidine_base_pattern = Chem.MolFromSmarts("n1[cH][cH][c,n][cH][cH]1")
    if not mol.HasSubstructMatch(pyrimidine_base_pattern):
       return False, "No pyrimidine base detected"

    # Deoxyribose sugar with a C-C-C-C-O ring system, where C2 has no OH.
    # We are using [C] instead of [C@H] or [C@@H] because the stereochemistry can vary, 
    # we are using the [] to indicate generic carbon atom.
    deoxyribose_pattern = Chem.MolFromSmarts("[C]1[C]([C]([C](O)CO)[C])O1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
         return False, "No deoxyribose sugar detected"


    # Check if there is ribose sugar instead of deoxyribose - if this pattern matches, it's not deoxy
    ribose_pattern = Chem.MolFromSmarts("[C]1[C](O[C]([C](O)CO)[C])O1")
    if mol.HasSubstructMatch(ribose_pattern):
        return False, "Ribose detected, not deoxyribose"


    #Verify glycosidic bond - N of pyrimidine is linked to C1' of sugar
    # This pattern looks for a nitrogen connected to a carbon, that is part of the ribose structure
    glycosidic_bond_pattern = Chem.MolFromSmarts("n~[C]1[C]([C]([C](O)CO)[C])O1")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond detected between pyrimidine and sugar"
        
    return True, "Pyrimidine deoxyribonucleoside detected"