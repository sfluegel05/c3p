"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a chlorobenzene carrying four chloro groups
    at unspecified positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all aromatic rings
    for ring in mol.GetRingInfo().AtomRings():
        # Ensure it's an aromatic ring
        if all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring):
            # Count chlorine atoms in the ring
            chloro_count = sum(1 for atom_idx in ring 
                               if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 17)
            if chloro_count == 4:
                return True, f"Contains an aromatic ring with four chlorine atoms"
    
    return False, "Does not match the tetrachlorobenzene pattern"