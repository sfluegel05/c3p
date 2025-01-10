"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound has bonds between metalloid atoms (e.g., arsenic) and carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find arsenic atoms
    arsenic_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 33]
    if not arsenic_atoms:
        return False, "No metalloid atoms (e.g., arsenic) found"
    
    # Check for bonds between arsenic and carbon
    for atom_idx in arsenic_atoms:
        for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                return True, "Contains metalloid-carbon bond: an organometalloidal compound"
    
    return False, "No metalloid-carbon bonds found"

# Example usage:
# smiles = "C[As](O)(O)=O"
# result, reason = is_organometalloidal_compound(smiles)
# print(result, reason)  # Expected: True, "Contains metalloid-carbon bond: an organometalloidal compound"