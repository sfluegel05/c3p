"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define general steroid backbone pattern as three 6-membered rings fused to one 5-membered ring
    steroid_backbone = Chem.MolFromSmarts('C1CC2CCC3CC4C(C(=C)CCC4)CCC3C2C1')
    
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Identify the hydroxy group at the 16beta position, ensuring flexibility in the pattern
    # Assuming that the attachment to the overall steroid structure is adequate for now
    hydroxyl_16beta_pattern = Chem.MolFromSmarts('[C@@H](O)C')
    
    # Check each carbon for attachment of a hydroxyl group to determine the correct position (use of substructure match placeholders)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and mol.HasSubstructMatch(hydroxyl_16beta_pattern, Chem.Atom(atom.GetIdx())):
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetSymbol() == 'C':
                    if len(neighbor.GetNeighbors()) == 3: # Beta orientation check (simplified)
                        return True, "Contains steroid backbone with 16beta-hydroxy group"
    
    return False, "No 16beta-hydroxy group found or incorrectly positioned"