"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    This class contains a glycerol backbone with one acyl group (ester bond) attached,
    while the other hydroxyl groups can be free or substituted.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Recognize glycerol backbone: three carbon chain with three hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2](O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Check for ester linkage with a long carbon chain (acyl group)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CH2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, "Need exactly one ester linkage for a monoacylglycerol"

    # Verify molecular structure according to typical monoacylglycerol structure
    # Count carbons to ensure there is a valid acyl chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if carbon_count < 10:  # Simple check for acyl chain length
        return False, f"Too few carbons ({carbon_count}) for a typical acyl chain, expected a longer chain"
    
    # Success
    return True, "Contains glycerol backbone with one acyl group attached via ester bond"