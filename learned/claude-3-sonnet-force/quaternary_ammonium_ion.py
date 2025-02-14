"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:53054 quaternary ammonium ion
A derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen
have been replaced with univalent (usually organyl) groups.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdqueries

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find positively charged nitrogen atoms with 4 substituents
    quat_n_query = Chem.MolFromSmarts("[N+;H0;D4]")
    quat_n_candidates = mol.GetSubstructMatches(quat_n_query)
    
    for candidate_idx in quat_n_candidates:
        candidate_atom = mol.GetAtomWithIdx(candidate_idx)
        
        # Check if all substituents are univalent (usually organyl) groups
        substituents = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in candidate_atom.GetNeighbors()]
        if all(is_univalent_group(sub_atom, mol) for sub_atom in substituents):
            return True, "Contains a positively charged quaternary nitrogen with 4 univalent substituents"
    
    return False, "No quaternary ammonium ion found"

def is_univalent_group(atom, mol):
    """
    Checks if the given atom is part of a univalent (usually organyl) group.

    Args:
        atom (Atom): The atom to check.
        mol (Mol): The molecule containing the atom.

    Returns:
        bool: True if the atom is part of a univalent group, False otherwise.
    """
    # Check if the atom is part of an alkyl, aryl, or other organic group
    if atom.IsInRingSize(5) or atom.IsInRingSize(6):
        return True  # Aryl group
    
    neighbors = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in atom.GetNeighbors()]
    if all(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):
        return True  # Alkyl group
    
    # Add additional checks for other univalent groups if needed
    
    return False