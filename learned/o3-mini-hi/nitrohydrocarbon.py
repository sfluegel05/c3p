"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: Nitrohydrocarbon
Definition: A C-nitro compound that is a hydrocarbon in which one or more of the hydrogens 
has been replaced by nitro groups.
"""

from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon (with only C as the backbone atoms) that has one or more  
    nitro ([N+](=O)[O-]) substituents attached to carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string to get a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a nitro group.
    # This pattern looks for an N with a positive charge bonded to two oxygens, one with a negative charge.
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    
    if not nitro_matches:
        return False, "No nitro groups found"
    
    # Check that at least one nitro group is attached to a carbon atom.
    # In a nitro group, the pattern [N+](=O)[O-] does not include the carbon atom that the nitro group
    # is attached to. Thus, for each match, look at the N atom (first in the match) and examine its neighbors.
    valid_nitro_found = False
    for match in nitro_matches:
        nitro_n_idx = match[0]  # index of the nitrogen in the nitro group
        nitro_n = mol.GetAtomWithIdx(nitro_n_idx)
        # Check all neighbors of the nitro nitrogen.
        # Two of them will be oxygens as in the pattern, the third (if present) should be the carbon of the hydrocarbon.
        for neighbor in nitro_n.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # carbon
                valid_nitro_found = True
                break
        if valid_nitro_found:
            break
    
    if not valid_nitro_found:
        return False, "No nitro group found attached to a carbon atom"
    
    # Next, gather all atom indices that are part of any nitro group.
    # We allow nitro group atoms (N and O obviously) to be present,
    # but the remainder of the molecule (the hydrocarbon scaffold) should be carbon atoms.
    nitro_atom_indices = set()
    for match in nitro_matches:
        for idx in match:
            nitro_atom_indices.add(idx)
    
    # Now check every heavy atom that is NOT part of a nitro group.
    for atom in mol.GetAtoms():
        # Only check heavy atoms (ignore implicit hydrogens)
        if atom.GetIdx() in nitro_atom_indices:
            continue
        # For the rest of the molecule, the atom should be carbon (atomic number 6).
        if atom.GetAtomicNum() != 6:
            return False, f"Found non-hydrocarbon atom ({atom.GetSymbol()}) in the carbon scaffold"
    
    # If all tests pass, we classify the molecule as a nitrohydrocarbon.
    return True, "Molecule is a nitrohydrocarbon: contains nitro groups attached to an otherwise pure carbon scaffold"
    
# Example usage:
if __name__ == "__main__":
    # Try a couple of examples:
    examples = {
        "3,7-Dinitrofluoranthene": "[O-][N+](=O)c1cccc-2c1-c1cccc3c(ccc-2c13)[N+]([O-])=O",
        "1-nitroheptane": "[O-][N+](=O)CCCCCCC",
        "Non-nitro hydrocarbon": "CCCCCC"
    }
    for name, smiles in examples.items():
        result, reason = is_nitrohydrocarbon(smiles)
        print(f"{name}: {result}. Reason: {reason}")