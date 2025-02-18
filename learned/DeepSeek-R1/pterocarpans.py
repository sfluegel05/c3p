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
    
    # Updated SMARTS pattern for benzofuro[3,2-c]chromene core without stereochemistry
    # Matches two benzene rings fused via oxygen-containing bridges
    pterocarpan_core = Chem.MolFromSmarts(
        "[O]1C2C3=C(C1)C=CC=C3COC4=C2C=CC=C4"
    )
    
    # Alternative pattern accounting for aromatic bonds and possible substitutions
    alt_core = Chem.MolFromSmarts(
        "c12c(ccc3c1OCC1c3cccc1)CO2"
    )
    
    # Check for core structure presence
    if mol.HasSubstructMatch(pterocarpan_core) or mol.HasSubstructMatch(alt_core):
        return True, "Contains benzofurochromene skeleton"
    
    # Additional check for fused ring system with oxygen bridges
    # This pattern looks for the fused benzofuro[3,2-c]chromene system
    fused_core = Chem.MolFromSmarts(
        "[O]1C2C3=C1C=CC=C3C4=C2COC=C4"
    )
    if mol.HasSubstructMatch(fused_core):
        return True, "Contains fused benzofurochromene structure"
    
    return False, "Lacks pterocarpan core structure"