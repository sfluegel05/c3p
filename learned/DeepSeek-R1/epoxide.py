"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:52092 epoxide
"""
from rdkit import Chem
from rdkit.Chem import GetSSSR

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is a cyclic ether with a 3-membered ring containing an oxygen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings in the molecule
    rings = mol.GetRingInfo().AtomRings()
    
    # Check each ring for 3-membered size and presence of oxygen
    for ring in rings:
        if len(ring) == 3:
            # Check if any atom in the ring is oxygen
            oxygen_in_ring = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)
            if oxygen_in_ring:
                return True, "Contains a 3-membered ring with oxygen (epoxide)"
    
    return False, "No 3-membered ring with oxygen found"