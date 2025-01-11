"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for sulfonic acid group (S(=O)(=O)O)
    sulfonic_acid_pattern = Chem.MolFromSmarts("S(=O)(=O)O")
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group found"
    
    # Check for C-S bond specifically attached to sulfonic group
    sulfur_atom_matches = mol.GetSubstructMatches(sulfonic_acid_pattern)
    for match in sulfur_atom_matches:
        sulfur_atom_idx = match[0]  # Sulfur atom index
        sulfur_atom = mol.GetAtomWithIdx(sulfur_atom_idx)
        # Find carbon atoms bonded to this sulfur
        c_s_bonded = any(bonded_atom.GetAtomicNum() == 6 for bonded_atom in sulfur_atom.GetNeighbors())
        if c_s_bonded:
            break
    else:
        return False, "No carbon-sulfur bond found involving sulfonic acid"
    
    # Check for lipid-like features (long carbon chains)
    # This is a simplified check and could be expanded for more precise lipid structure
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:  # Assuming lipid has at least 20 carbons
        return False, "Too few carbons to suggest a lipid-like structure"

    return True, "Contains sulfonic acid group bonded to a carbon in a lipid structure"