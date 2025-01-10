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

    # Extend isoprene pattern to include cis/trans configurations detected via varied SMILES junction
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=C[CH]=C[CH2]")  # General pattern for isoprene units
    
    # Special check for polyprenol chains
    long_isoprene_chain = Chem.MolFromSmarts("(C=C\C=C\C=C)+")

    # Find isoprene matches
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    long_chain_matches = mol.HasSubstructMatch(long_isoprene_chain)
    
    # Check if there is a sensible amount linked as a single chain for polyprenol property
    if len(isoprene_matches) < 4 or not long_chain_matches:
        return False, "Insufficient linked isoprene units for a polyprenol chain (need at least 4 in sequence)"

    # Recognize phosphate groups;
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])OP(=O)(O)O")
    
    # Find phosphate matches
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)

    if not phosphate_matches and not diphosphate_matches:
        return False, "No phosphate or diphosphate group found"
    
    # Ensure the phosphate or diphosphate group is linked at what can be considered the terminal of the isoprene chain
    terminal_atoms = {match[0] for match in isoprene_matches}.union({match[-1] for match in isoprene_matches})

    # Test for correct phosphate linkage
    attached_phosphate = False
    
    for match in phosphate_matches + diphosphate_matches:
        phosphate_end = {bond.GetBeginAtomIdx() for bond in mol.GetBonds() if bond.IsInRingSize(5)}
        if any(atom in terminal_atoms for atom in match) or any(atom in phosphate_end for atom in match):
            attached_phosphate = True
            break

    if attached_phosphate:
        return True, "Valid polyprenol phosphate with correct phosphate attachment to isoprene chain"
    
    return False, "Phosphate group is incorrectly positioned in the isoprene chain"