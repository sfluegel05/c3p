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
    This corresponds to a chromen-4-one core with an aryl group at position 3.

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
    
    # Correct SMARTS pattern for chromen-4-one core with atom mapping
    chromen4one_smarts = """
        [cH]1-[cH]-[c]2-[o]-[c](=[O])-[cH]-[c]2-[cH]-[cH]1
        """
    chromen4one_mol = Chem.MolFromSmarts(chromen4one_smarts)
    if chromen4one_mol is None:
        return False, "Invalid SMARTS pattern for chromen-4-one core"
    
    # Find substructure matches of the chromen-4-one core
    matches = mol.GetSubstructMatches(chromen4one_mol)
    if not matches:
        return False, "No chromen-4-one core found"
    
    # For each match, check for an aryl substituent at position 3
    for match in matches:
        # Atom indices of the matched core in the molecule
        core_atom_indices = list(match)
        
        # Create a map from pattern atom indices to molecule atom indices
        # This depends on the order of atoms in the SMARTS pattern
        # The position of atom indices in the SMARTS pattern corresponds to the atom indices in the match
        # We need to identify position 3 in the core
        # Based on the SMARTS pattern, position 3 corresponds to the 6th atom (index 5)
        pos3_atom_idx = match[5]  # 0-based indexing
        
        # Get the atom at position 3 in the molecule
        pos3_atom = mol.GetAtomWithIdx(pos3_atom_idx)
        
        # Check neighbors of position 3 atom
        neighbors = pos3_atom.GetNeighbors()
        for nbr in neighbors:
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in core_atom_indices:
                # Check if neighbor atom is part of an aromatic ring not in the core
                # Use GetRingInfo()
                rings = mol.GetRingInfo().AtomRings()
                is_aromatic_ring = False
                for ring in rings:
                    if nbr_idx in ring and not set(ring).issubset(core_atom_indices):
                        # Check if all atoms in the ring are aromatic
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