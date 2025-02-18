"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: CHEBI:18214 isoflavones
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone has a 3-aryl-4H-chromen-4-one (3-aryl-4H-chromen-4-one) skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core structure pattern: 3-aryl-4H-chromen-4-one
    # Chromen-4-one is a benzene fused to pyrone (O=C1C=CC2=C1C=CC=C2)
    # The aryl group is attached to the 3-position (third carbon in pyrone)
    core_pattern = Chem.MolFromSmarts("c1ccc2c(c1)occc(=O)c2-[c]")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing 3-aryl-4H-chromen-4-one core"

    # Verify the aryl group is part of an aromatic ring (at least 6-membered)
    matches = mol.GetSubstructMatches(core_pattern)
    for match in matches:
        # The last atom in the match is the connecting carbon of the aryl group
        aryl_start = match[-1]
        # Check if the substituent is part of an aromatic ring
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if aryl_start in ring:
                # Check if the ring is aromatic and has at least 6 atoms
                if len(ring) >= 6 and all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
                    return True, "Contains 3-aryl-4H-chromen-4-one skeleton"
        
    return False, "Aryl substituent is not part of an aromatic ring"