"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: CHEBI:????? corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid contains a corrin nucleus: four pyrrole-like rings connected in a macrocycle
    with three =C- groups and one direct C-C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the corrin nucleus SMARTS pattern
    corrin_smarts = Chem.MolFromSmarts("C1C-C=C-C=N-C(-CC)-C=C-N-C(-CC)-C=C-N-C(-CC)-C1-N1")
    
    # Check if the corrin nucleus is present
    if mol.HasSubstructMatch(corrin_smarts):
        return True, "Contains corrin nucleus with four pyrrole rings connected by three =C- groups and one C-C bond"
    
    return False, "No corrin nucleus detected"