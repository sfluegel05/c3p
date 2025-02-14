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

    # Aniline feature: Phenyl group attached to a nitrogen directly forming an amide
    aniline_pattern = Chem.MolFromSmarts("c1ccc(cc1)NC(=O)")
    # Allows variation in acyl group by expanding the pattern to include any carbonyl attached to N
    extended_pattern = Chem.MolFromSmarts("c1ccc(cc1)N=[OD1]")

    # Check if the molecule contains an aniline-like substructure
    if mol.HasSubstructMatch(aniline_pattern):
        return True, "Contains aniline-like aromatic amide structure"
    elif mol.HasSubstructMatch(extended_pattern):
        return True, "Contains extended aniline-like aromatic amide structure"
    else:
        return False, "Does not match anilide structural criteria"