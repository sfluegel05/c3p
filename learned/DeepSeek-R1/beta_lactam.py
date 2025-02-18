"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam contains a four-membered ring with an amide bond (N-C=O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for four-membered ring with amide bond (N-C=O)
    # [nH0] matches nitrogen in a ring with no hydrogen (common in fused rings)
    # The pattern looks for a 4-membered ring where N is adjacent to C=O
    beta_lactam_pattern = Chem.MolFromSmarts("[nH0]1C(=O)CC1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No four-membered lactam ring found"

    # Verify ring size - ensure matched atoms form exactly 4-membered ring
    matches = mol.GetSubstructMatches(beta_lactam_pattern)
    for match in matches:
        ring_atoms = list(match)
        ring = Chem.GetRingInfo(mol).AtomRings()
        for r in ring:
            if set(ring_atoms).issubset(r) and len(r) == 4:
                return True, "Contains four-membered lactam ring with amide bond"

    return False, "No valid four-membered lactam ring detected"