"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: 1-monoglyceride (CHEBI:17408)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride has a glycerol backbone with a single acyl group at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 1-monoglyceride core structure:
    # Glycerol backbone with ester at position 1 and hydroxyls at 2 and 3
    # Explicit connectivity to avoid partial matches:
    # [CH2]-[CH]-[CH2] where first CH2 has ester, middle CH has OH, last CH2 has OH
    glyceride_pattern = Chem.MolFromSmarts(
        "[CH2]([OX2][C](=O)*)-[CH]([OH])-[CH2]([OH])"
    )
    alt_pattern = Chem.MolFromSmarts(
        "[CH2]([OH])-[CH]([OH])-[CH2]([OX2][C](=O)*)"
    )

    # Check both possible orientations (SMILES direction matters)
    matches = mol.GetSubstructMatches(glyceride_pattern)
    alt_matches = mol.GetSubstructMatches(alt_pattern)
    
    if not (matches or alt_matches):
        return False, "Glycerol backbone with ester at position 1 not found"

    # Verify there's exactly one ester group in the entire molecule
    ester_pattern = Chem.MolFromSmarts("[OX2][C](=O)")
    ester_count = len(mol.GetSubstructMatches(ester_pattern))
    if ester_count != 1:
        return False, f"Found {ester_count} ester groups (must be exactly 1)"

    # Ensure no other substituents on glycerol backbone
    # Check that all oxygen atoms are either in hydroxyls or the single ester
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            if not (atom.GetDegree() == 2 and atom.GetTotalNumHs() == 0):  # Ester oxygen
                if not (atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1):  # Hydroxyl
                    return False, "Unexpected oxygen substituents"

    # Check molecular formula constraints (C3H8O3 + acyl group)
    # Total oxygens should be 3 (2 OH + 1 ester oxygen)
    o_count = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    if o_count != 3 + 1:  # Glycerol (3O) + ester adds 1O from the acyl
        return False, f"Unexpected oxygen count: {o_count}"

    return True, "1-monoglyceride with acyl group at position 1"