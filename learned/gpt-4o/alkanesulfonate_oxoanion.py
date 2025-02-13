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

    # Look for alkanesulfonate (-SO3^-) structure: S([O-])(=O)=O
    sulfonate_pattern = Chem.MolFromSmarts("S([O-])(=O)=O")
    matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not matches:
        return False, "No sulfonate group S([O-])(=O)=O found"

    # Check if the sulfonate group is attached specifically to a saturated carbon
    # (i.e., in an alkane context) with permissible R groups (hydrogen, carbon, others)
    for match in matches:
        sulfur_atom = mol.GetAtomWithIdx(match[0])
        
        # Check for carbon neighbors of the sulfur
        carbon_neighbors = [n for n in sulfur_atom.GetNeighbors() if n.GetAtomicNum() == 6]
        for carbon in carbon_neighbors:
            # Verify if the carbon is part of an alkane (saturated with hydrogens or other single-bonded groups)
            # A crude check is to ensure the carbon has 4 bonds, typically indicating saturation
            if len(carbon.GetNeighbors()) == 4:
                return True, "Contains a sulfonate group attached to a saturated alkane carbon or similar R group"

    return False, "Sulfonate group not attached to a saturated alkane carbon or appropriate R group"