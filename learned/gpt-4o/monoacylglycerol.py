"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol comprises a glycerol backbone with one acyl group
    (ester bond) and potentially other substituents on the two remaining hydroxyls.

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

    # Recognize glycerol-like backbone: three carbon chain with potential substitutes
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4](O)[CH2X4](O) |$(*,A1)&&(*,A2)|")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol-like backbone found"
    
    # Check for ester linkage with only one acyl group attached
    ester_pattern = Chem.MolFromSmarts("O=C(O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 1 for a monoacylglycerol"

    # Ensure verified ester is with a long chain typical for acyl groups
    for atom in mol.GetAtoms():
        if atom.GetIdx() == ester_matches[0][2] and atom.GetDegree() > 1:
            if sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6) < 10:
                return False, "Acyl group carbon chain seems too short"

    return True, "Contains glycerol backbone with one acyl group attached via ester bond"