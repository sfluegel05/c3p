"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:18357 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine contains a catechol moiety (benzene ring with adjacent hydroxyl groups)
    and an aminoethyl side chain attached to the benzene ring, including derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define catechol moiety pattern (benzene ring with two adjacent hydroxyl groups)
    catechol_pattern = Chem.MolFromSmarts('a1aa(O)aa(O)aa1')
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol moiety found (benzene ring with adjacent hydroxyl groups)"
    
    # Define aminoethyl side chain pattern attached to aromatic ring
    aminoethyl_pattern = Chem.MolFromSmarts('[NX3][CX4][CX4H2][c]')
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl side chain attached to aromatic ring found"
    
    # Ensure that the catechol moiety and aminoethyl side chain are connected
    # Get matches for catechol moiety
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    # Get matches for aminoethyl side chain
    aminoethyl_matches = mol.GetSubstructMatches(aminoethyl_pattern)
    
    # Check connectivity
    for cat_match in catechol_matches:
        catechol_atoms = set(cat_match)
        for amine_match in aminoethyl_matches:
            if any(neighbor.GetIdx() in catechol_atoms for atom_idx in amine_match for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors()):
                return True, "Contains catechol moiety with aminoethyl side chain attached"
    
    return False, "Catechol moiety and aminoethyl side chain not properly connected"