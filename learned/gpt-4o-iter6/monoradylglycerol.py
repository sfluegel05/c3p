"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol with one acyl, alkyl, or alk-1-enyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for glycerol backbone
    # [O][C][C][O] and [O][C][C][O] are interchangeable to allow flexibility in attachment points.
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for substituents (acyloxy, alkoxy, alkenyl) attached to glycerol
    substituent_patterns = [
        Chem.MolFromSmarts("C(=O)OCC"),  # Acyl group attached to an ether
        Chem.MolFromSmarts("C=C-OCC"),  # Alkenyl group via ether
        Chem.MolFromSmarts("COCC")      # Additional ethers or alkyl-type ranges
    ]

    # Count non-glycerol chain attachments
    substituent_count = sum(mol.HasSubstructMatch(pattern) for pattern in substituent_patterns)

    # Exactly one large substituent characterizes monoradylglycerol
    if substituent_count != 1:
        return False, f"Expected 1 substituent group, found {substituent_count}"

    return True, "Contains glycerol backbone with one acyl, alkyl, or alk-1-enyl substituent"