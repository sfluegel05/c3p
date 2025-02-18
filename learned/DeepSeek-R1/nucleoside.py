"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside consists of a ribose/deoxyribose sugar connected via an N-glycosidic bond to a nucleobase.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of a sugar moiety (ribose/deoxyribose pattern)
    # Generic furanose pattern with possible deoxy at C2 and hydroxymethyl (C5)
    sugar_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H](O[C@@H]1CO)O)O)")
    deoxy_sugar_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H](O[C@@H]1CO)O)C)")
    if not (mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(deoxy_sugar_pattern)):
        return False, "No ribose/deoxyribose sugar detected"

    # Check for N-glycosidic bond (sugar C1 connected to nucleobase N)
    glycosidic_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@H](O[C@@H]1CO)O)-[N]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No N-glycosidic bond detected"

    # Check for nucleobase (purine or pyrimidine derivatives)
    purine_core = Chem.MolFromSmarts("n1c2nc[nH]c2nc1")  # Purine core (adenine/guanine)
    pyrimidine_core = Chem.MolFromSmarts("n1c(=O)[nH]c(=O)cc1")  # Pyrimidine core (uracil)
    pyrimidine_alt = Chem.MolFromSmarts("n1c(=O)ccnc1")  # Cytosine/thymine variants
    
    if (mol.HasSubstructMatch(purine_core) or 
        mol.HasSubstructMatch(pyrimidine_core) or 
        mol.HasSubstructMatch(pyrimidine_alt)):
        return True, "Contains sugar and nucleobase with N-glycosidic bond"
    
    # Handle modified bases with different substitution patterns
    modified_purine = Chem.MolFromSmarts("[n]1[c][n][c]2[c][n][c][c]12")  # Modified purine scaffold
    modified_pyrimidine = Chem.MolFromSmarts("[n]1[c][c][c][c][c]1=O")    # Modified pyrimidine
    
    if (mol.HasSubstructMatch(modified_purine) or 
        mol.HasSubstructMatch(modified_pyrimidine)):
        return True, "Contains modified nucleobase structure"
    
    return False, "No valid nucleobase detected"