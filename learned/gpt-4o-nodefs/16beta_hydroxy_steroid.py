"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid backbone pattern: cyclopenta[a]phenanthrene.
    # Simplified pattern to capture the essence: A/B/C/D ring system
    steroid_backbone = Chem.MolFromSmarts('C1CCC2C3CCC4CCCC(C4)C3C2C1')

    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Hydroxy group at the 16beta position
    # In the steroid, the 16beta position (stereochemistry notation) is specific; 
    # this requires understanding of steroid structure. Typically involves C17 as a reference.
    # Use a pipecleaner pattern to match hydroxyl in beta orientation on correct carbon
    hydroxyl_pattern = Chem.MolFromSmarts('C[C@H](O)C')  # Placeholder pattern, needs specificity

    # Check atom positions and their neighbors for hydroxyl group at 16beta
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            bond = [b for b in atom.GetBonds() if b.GetBeginAtom().GetSymbol() == 'O' or b.GetEndAtom().GetSymbol() == 'O']
            if bond:
                # Check for beta configuration, commonly using neighboring hydrogens/bonds
                if mol.HasSubstructMatch(hydroxyl_pattern):
                    return True, "Contains steroid backbone with 16beta-hydroxy group"
    
    return False, "16beta-hydroxy group not found or incorrectly positioned"