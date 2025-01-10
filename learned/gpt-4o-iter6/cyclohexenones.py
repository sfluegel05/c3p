"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as any six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expanded SMARTS to account for variations in cyclohexenone structures
    # 1. It caters for any carbon at different positions and double bonds distribution
    # 2. Allows substituents and is flexible on position of the ketone
    # 3. *C indicates any atom can be attached to either of the carbons
    cyclohexenone_pattern = Chem.MolFromSmarts("C1=CC(=O)[C,C][C,C][C,C]1")
    
    # Check for the cyclohexenone pattern
    if not mol.HasSubstructMatch(cyclohexenone_pattern):
        return False, "No cyclohexenone structure found"
    
    # Check the number of carbon atoms in the ring to ensure it remains six-membered
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if len(ring) == 6 and any(atom.GetHybridization() == Chem.HybridizationType.SP2 for atom in ring_atoms):
            return True, "Cyclohexenone structure identified with appropriate ring and ketone"
    
    return False, "No six-membered alicyclic ketone with one double bond found in the ring"