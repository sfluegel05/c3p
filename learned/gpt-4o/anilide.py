"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is defined as any aromatic amide obtained by acylation of aniline.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Anilide feature: Phenyl group directly attached to a nitrogen forming an amide (expanded for substitution on the phenyl)
    anilide_patterns = [
        Chem.MolFromSmarts("c1ccccc1NC(=O)"),  # Unsubstituted phenylamide
        Chem.MolFromSmarts("c1cc(c(cc1)*)NC(=O)")  # Substituted phenyl containing amide
    ]
    
    # Check for any pattern match that will qualify it as an anilide
    for pattern in anilide_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains an aniline-like aromatic amide structure"

    return False, "Does not match anilide structural criteria"