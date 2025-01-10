"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside is an N-glycosyl compound that has both a nucleobase and a sugar moiety (ribose or deoxyribose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nucleobase patterns (adenine, guanine, thymine, cytosine, uracil, xanthine, and their modifications)
    nucleobase_patterns = [
        Chem.MolFromSmarts("[nH]1cnc2c1nc[nH]2"),  # Adenine
        Chem.MolFromSmarts("[nH]1c(=O)[nH]c2c1ncn2"),  # Guanine
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1cc[nH]2"),  # Cytosine
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1ccc[nH]2"),  # Uracil
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1c(C)cc[nH]2"),  # Thymine
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1c(O)cc[nH]2"),  # Xanthine
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1c(N)cc[nH]2"),  # Modified Cytosine
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1c(F)cc[nH]2"),  # Modified Uracil
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1c(Br)cc[nH]2"),  # Modified Uracil
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1c(Cl)cc[nH]2"),  # Modified Uracil
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1c(CO)cc[nH]2"),  # Modified Uracil
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1c(OC)cc[nH]2"),  # Modified Uracil
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1c(SC)cc[nH]2"),  # Modified Uracil
        Chem.MolFromSmarts("[nH]1c(=O)nc2c1c(NC)cc[nH]2")   # Modified Uracil
    ]

    # Check for presence of a nucleobase
    has_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not has_nucleobase:
        return False, "No nucleobase found"

    # Define sugar patterns (ribose or deoxyribose, including modifications)
    ribose_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@H](O)[C@H]1O")  # Ribose
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@H](O)[C@H]1")  # Deoxyribose
    modified_ribose_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@H](O)[C@H]1[O,N,S]")  # Modified Ribose
    modified_deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@H](O)[C@H]1[O,N,S]")  # Modified Deoxyribose

    # Check for presence of a sugar moiety
    has_sugar = (mol.HasSubstructMatch(ribose_pattern) or 
                 mol.HasSubstructMatch(deoxyribose_pattern) or 
                 mol.HasSubstructMatch(modified_ribose_pattern) or 
                 mol.HasSubstructMatch(modified_deoxyribose_pattern))
    if not has_sugar:
        return False, "No sugar moiety (ribose or deoxyribose) found"

    # Check for N-glycosidic bond between nucleobase and sugar
    glycosidic_bond_pattern = Chem.MolFromSmarts("[NX3][C@H]1O[C@H](CO)[C@H](O)[C@H]1[O,N,S]")  # N-glycosidic bond
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No N-glycosidic bond found between nucleobase and sugar"

    return True, "Contains a nucleobase and a sugar moiety connected via an N-glycosidic bond"