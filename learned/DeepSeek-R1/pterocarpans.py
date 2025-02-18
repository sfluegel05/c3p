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
    
    # Core pattern matching the benzofuro[3,2-c]chromene skeleton
    # Accounts for fused benzene rings with oxygen bridge and dihydro arrangement
    core_pattern = Chem.MolFromSmarts(
        "[#6]12[#6](:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#8]-[#6]-[#6]-1)[#8][#6](:[#6]:[#6]:[#6]:[#6]:[#6]:21)"
    )
    
    # Alternative pattern using aromatic notation for better matching
    alt_core = Chem.MolFromSmarts(
        "c12c3c(ccc4c3OCCc3c4cccc3)COc1cccc2"
    )
    
    # Check for core structure presence
    if mol.HasSubstructMatch(core_pattern) or mol.HasSubstructMatch(alt_core):
        return True, "Contains benzofurochromene skeleton"
    
    # Final check with relaxed pattern focusing on oxygen bridges
    relaxed_core = Chem.MolFromSmarts(
        "c1ccc2c3c1OCCc1c3cccc1CO2"
    )
    if mol.HasSubstructMatch(relaxed_core):
        return True, "Contains fused benzofurochromene structure"
    
    return False, "Lacks pterocarpan core structure"