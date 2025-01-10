"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    Wax esters consist of long-chain fatty acids and fatty alcohols connected through ester linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for ester linkage pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if not ester_matches:
        return False, "No ester linkage found"
    
    for match in ester_matches:
        c_atom, o1_atom, o2_atom = match

        # Check if each side of the ester linkage has at least 14 carbons
        left_chain_length = get_chain_length(mol, c_atom, {o1_atom, o2_atom})
        right_chain_length = get_chain_length(mol, o2_atom, {c_atom})

        # Consider only the side from the ester C flipped direction (RCO-O)
        if left_chain_length >= 14 and right_chain_length >= 14:
            # Further check for additional elements that could wrongly contribute
            # Exclude if molecule contains complex groups indicative of non-wax esters
            phosphate_group = Chem.MolFromSmarts("P(=O)(O)O")
            if mol.HasSubstructMatch(phosphate_group):
                continue
            
            amine_group = Chem.MolFromSmarts("N")
            if mol.HasSubstructMatch(amine_group) and is_seed - complex_lipid_like(mol):  # Function checking for complex lipid
                continue

            return True, "Contains ester linkage with sufficient long carbon chains, matching wax ester"

    return False, "Ester linkages found, but chains are not sufficiently long or match complex non-wax ester"

def get_chain_length(mol, start_atom_idx, visited):
    """Extend chain length calculation from a given atom, skipping visited atoms."""
    length = 0
    to_visit = {start_atom_idx}
    
    while to_visit:
        current_atom_idx = to_visit.pop()
        current_atom = mol.GetAtomWithIdx(current_atom_idx)
        if current_atom.GetAtomicNum() == 6:  # Count only carbons
            length += 1
        
        for nbr in current_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in visited:
                visited.add(nbr_idx)
                to_visit.add(nbr_idx)
    
    return length

def is_complex_lipid_like(mol):
    """Placeholder to further check for complex molecules like phospholipids."""
    # Implement checks for additional structural complexity here, if needed
    return False