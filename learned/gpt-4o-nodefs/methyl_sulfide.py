"""
Classifies: CHEBI:86315 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is characterized by a direct bond between a terminal methyl group
    and a sulfur atom (CH3-S-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a terminal methylthio group pattern (CH3-S and no other group on the sulfur)
    # This pattern specifies that the S atom connects to a CH3 group,
    # and considers the context to avoid branching structures at the sulfur site.
    methylthio_strict_pattern = Chem.MolFromSmarts("[C]([H3])[S]")
    sulfur_pattern = Chem.MolFromSmarts("[S]")
    
    # Check if there is any methyl group attached as a terminal to a sulfur
    strict_matches = mol.HasSubstructMatch(methylthio_strict_pattern)
    
    # Ensure sulfur itself is available with only the methyl group attached
    if strict_matches and all(len(atom.GetNeighbors()) <= 2 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16):
        return True, "Contains methylthio group (CH3-S-) without branching"
    
    return False, "No suitable methylthio group found"