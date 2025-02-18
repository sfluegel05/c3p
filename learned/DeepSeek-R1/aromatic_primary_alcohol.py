"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:134307 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a hydroxyl group (-OH) attached to a primary carbon
    that is directly bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for primary alcohol attached to any aromatic atom
    # [a] is any aromatic atom, [CH2;D2] is methylene group with exactly two bonds (to OH and aromatic atom)
    pattern = Chem.MolFromSmarts('[a]-[CH2;D2]-[OH]')
    
    # Check for matches
    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No primary alcohol group adjacent to aromatic ring found"
    
    # Verify that the aromatic atom is part of a ring
    valid = False
    for match in matches:
        aromatic_atom_idx = match[0]
        atom = mol.GetAtomWithIdx(aromatic_atom_idx)
        if atom.IsInRing():
            valid = True
            break
    
    if valid:
        return True, "Primary alcohol group attached to aromatic carbon in a ring"
    else:
        return False, "Aromatic atom not part of a ring"