"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is likely classified as a semisynthetic derivative.
    It checks for the presence of natural skeletons with synthetic modifications.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely semisynthetic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Search for natural product core structures (e.g., steroids, lactones, etc.)
    steroid_pattern = Chem.MolFromSmarts("[#6]1-[#6]-[#6]2-[#6]-[#6]3-[#6]-[#6]4-[#6]-1-[#6]-[#6]2-[#6]-[#6]3-[#6]-[#6]4")
    lactone_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#8]:[#6]-1")
    
    # Check for synthetic modifications such as unusual bridging or extra functional groups
    unusual_fu_group = Chem.MolFromSmarts("[#6][C@](=[O,N])[O,N,c]")  # includes acetals, amide bridge
    protecting_group_pattern = Chem.MolFromSmarts("[#6][O][C]([#6])[#6]")  # typical protecting groups
    
    steroid_matches = mol.HasSubstructMatch(steroid_pattern)
    lactone_matches = mol.HasSubstructMatch(lactone_pattern)
    
    unusual_matches = mol.GetSubstructMatches(unusual_fu_group)
    protecting_group_matches = mol.GetSubstructMatches(protecting_group_pattern)
    
    # Evaluate matches and make a decision
    if (steroid_matches or lactone_matches) and (len(unusual_matches) > 0 or len(protecting_group_matches) > 0):
        reason = "Molecule has natural core with synthetic modifications: "
        reason += f"Unusual functional groups ({len(unusual_matches)}), " if len(unusual_matches) > 0 else ""
        reason += f"Protecting group patterns ({len(protecting_group_matches)})" if len(protecting_group_matches) > 0 else ""
        return True, reason.strip(", ")
    
    return False, "No characteristic semisynthetic derivative patterns identified"