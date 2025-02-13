"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect sulfonate group with [O-]: S([O-])(=O)=O or S([O-])(=O)(=O)
    sulfonate_pattern = Chem.MolFromSmarts("S([O-])(=O)(=O)")
    matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not matches:
        return False, "No sulfonate group S([O-])(=O)(=O) found"

    # Check connections of the sulfur atom to ensure it is bonded correctly
    for match in matches:
        sulfur_atom = mol.GetAtomWithIdx(match[0])
        
        # Check for carbon or permitted atom neighbors of sulfur; extending beyond just saturated carbon
        connected_atoms = [n for n in sulfur_atom.GetNeighbors() if n.GetAtomicNum() in [6, 1, 7, 8, 15]]
        if connected_atoms:
            return True, "Found an appropriately attached sulfonate group with a diverse R group"

    return False, "Sulfonate group not properly attached"