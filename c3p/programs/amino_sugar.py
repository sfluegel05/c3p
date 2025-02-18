"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is defined as any sugar having one or more alcoholic hydroxy groups
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Recognize sugar-like rings more broadly by considering other ring types
    sugar_patterns = [
        Chem.MolFromSmarts("C[C@H]1O[C@H](O)[C@@H](CO)[C@@H](N)C1"), # Generic 6-membered sugar ring with amino
        Chem.MolFromSmarts("C1OC([C@@H](O)C[C@H](N)O1)"), # 6-membered ring with an amine
        Chem.MolFromSmarts("C1OC(C[NH]C[C@@H](O)[C@@H]1)"), # 5-membered furanose-like with amine
    ]
    
    found_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(pattern):
            found_sugar = True
            break
    
    if not found_sugar:
        return False, "No typical sugar-like ring structure found"
    
    # Confirm presence of amino groups in potential substitution sites for hydroxyls
    nh_group = Chem.MolFromSmarts("[NH2,NH1,N#N]")
    if not mol.HasSubstructMatch(nh_group):
        return False, "No amino group replacing a hydroxy group found"
    
    return True, "Contains a sugar-like ring structure with one or more hydroxyl groups replaced by amino groups"