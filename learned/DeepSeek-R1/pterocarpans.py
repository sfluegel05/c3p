"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: CHEBI:26308 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans have a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # SMARTS pattern for the pterocarpan core skeleton:
    # Benzofuro[3,2-c]chromene with 6a,11a-dihydro (two adjacent CH2 groups in the chromene part)
    pterocarpan_core = Chem.MolFromSmarts(
        "[O]1c2ccccc2C2[C@@H]1C1c3c(cccc3)O[C@H]1C2"  # Accounts for stereochemistry in fused rings
    )
    
    # Check if the core structure is present
    if mol.HasSubstructMatch(pterocarpan_core):
        return True, "Contains benzofurochromene skeleton with 6a,11a-dihydro structure"
    
    # Alternative check without stereochemistry if the first fails
    alt_core = Chem.MolFromSmarts(
        "O1c2ccccc2C2C1C1c3c(cccc3)OC1C2"  # Non-stereospecific version
    )
    if mol.HasSubstructMatch(alt_core):
        return True, "Contains benzofurochromene skeleton (non-stereo)"
    
    return False, "Lacks pterocarpan core structure"