"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
"""
from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the chromen-4-one core (benzopyran-4-one)
    chromen4one_smarts = 'O=C1C=CC=CC2=CC=CC=C12'  # Chromen-4-one core
    chromen4one = Chem.MolFromSmarts(chromen4one_smarts)
    if chromen4one is None:
        return False, "Invalid SMARTS pattern for chromen-4-one core"
    
    # Find substructure matches of the chromen-4-one core
    matches = mol.GetSubstructMatches(chromen4one)
    if not matches:
        return False, "No chromen-4-one core found"
    
    # For each match, check for an aryl substituent at position 3
    for match in matches:
        # Atoms in chromen-4-one core
        core_atoms = set(match)
        
        # Map of atom indices in chromen4one_smarts to indices in the molecule
        atom_map = {i: idx for i, idx in enumerate(match)}
        
        # Atom at position 3 in chromen-4-one core (index 4 in chromen4one_smarts)
        pos3_idx_core = 4  # Indexing starts from 0 in SMARTS pattern
        pos3_atom_idx = atom_map.get(pos3_idx_core, None)
        if pos3_atom_idx is None:
            continue  # Cannot find position 3 atom, skip to next match
        pos3_atom = mol.GetAtomWithIdx(pos3_atom_idx)
        
        # Check neighbors of position 3 atom
        neighbors = pos3_atom.GetNeighbors()
        for nbr in neighbors:
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in core_atoms:
                # Check if neighbor atom is part of an aromatic ring not in core
                # Use GetRingInfo()
                rings = mol.GetRingInfo().AtomRings()
                is_aromatic_ring = False
                for ring in rings:
                    if nbr_idx in ring and not core_atoms.intersection(ring):
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                            is_aromatic_ring = True
                            break
                if is_aromatic_ring:
                    return True, "Contains isoflavonoid core with aryl substituent at position 3"
        # If no aryl substituent found for this match, continue to next match
    return False, "Isoflavonoid core found, but no aryl substituent at position 3"

# Example usage:
# smiles = 'Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O'  # Daidzein, an isoflavonoid
# result, reason = is_isoflavonoid(smiles)
# print(result, reason)