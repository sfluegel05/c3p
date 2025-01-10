"""
Classifies: CHEBI:30527 flavin
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    A flavin is a derivative of dimethylisoalloxazine with a substituent on the 10th position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavin, False otherwise
        str: Reason for classification
    """
    # Convert SMILES string to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible SMARTS pattern for the dimethylisoalloxazine core
    # Include options for bithiazole variability and connectivity issues
    core_pattern = Chem.MolFromSmarts("Cc1cnc2c(C)c3[nH]c(=O)n(c3=O)c2nc1C")
    
    # Check for core structure
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core flavin structure (dimethylisoalloxazine) missing"

    # Determine if there's a probable substituent at the nitrogen atom (common flavin 10th position)
    substituent_match = False
    core_matches = mol.GetSubstructMatches(core_pattern)
    
    for match in core_matches:
        # Check which pattern positions align with substitutable sites
        # Consider neighboring atoms excluding those accounted in the core signature
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == 'N':  # Target the Nitrogen of interest (assuming similar position: match[4])
                for neighbor in atom.GetNeighbors():
                    # Substituents are non-hydrogen atoms bonded to the nitrogen in question
                    if neighbor.GetAtomicNum() not in [1, 7]:  # 1 is hydrogen, 7 is other nitrogen parts
                        substituent_match = True
                        break
            if substituent_match:
                break            

    if not substituent_match:
        return False, "No substituent detected at the expected nitrogen position"
    
    return True, "Valid flavin structure detected"