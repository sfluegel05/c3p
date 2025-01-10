"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    A pterocarpan is characterized by a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the substructure pattern for the core of pterocarpans
    # SMARTS pattern for the 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene core
    pterocarpan_pattern = Chem.MolFromSmarts("O1[C@@]2([C@](C3=C1C=CC=4OCOC4C3)([H])[H])") 

    # Check if the core structure is present
    if not mol.HasSubstructMatch(pterocarpan_pattern):
        return False, "No pterocarpan core found"
    
    # Further checks can be added here for additional modifications common in pterocarpans,
    # such as specific hydroxylation patterns or attached isoflavonoid-like structures.
    
    return True, "Contains pterocarpan core structure"