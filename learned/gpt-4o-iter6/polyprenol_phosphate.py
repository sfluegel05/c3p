"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a polyprenol phosphate, False otherwise.
        str: Reason for classification or failure.
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern to capture both trans and cis-isoprene units: CH2=C-CH=CH2 and variations
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=C([CH3])[CH]=[CH2]")
    
    # Find isoprene matches
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 4:
        return False, "Insufficient isoprene units for a polyprenol chain (need at least 4)"

    # Recognize phosphate and diphosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)[O-]")
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])OP(=O)([O-])[O-]")
    
    # Find phosphate matches
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not phosphate_matches and not diphosphate_matches:
        return False, "No phosphate or diphosphate group found"
    
    # Ensure the phosphate or diphosphate group is at the end of the isoprene chain
    terminal_atoms = {match[0] for match in isoprene_matches}.union({match[-1] for match in isoprene_matches})
    attached_phosphate = False
    for match in phosphate_matches + diphosphate_matches:
        if any(atom in terminal_atoms for atom in match):
            attached_phosphate = True
            break

    if attached_phosphate:
        return True, "Valid polyprenol phosphate discovered"
    
    return False, "Phosphate group is not correctly positioned"