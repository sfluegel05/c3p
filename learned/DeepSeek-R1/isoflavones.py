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
    An isoflavone has a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton.

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
    # The aryl group is attached to the 3-position (carbon adjacent to oxygen in pyrone)
    core_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)cc2-[c]")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing 3-aryl-4H-chromen-4-one core"

    # Verify the aryl group is aromatic (at least a benzene ring)
    # Get the matched atoms to check the aryl substituent
    matches = mol.GetSubstructMatches(core_pattern)
    for match in matches:
        # The last atom in the match is the connecting carbon of the aryl group
        aryl_start = match[-1]
        # Check if the aryl group is part of an aromatic ring
        aryl_atoms = [aryl_start]
        # Traverse the connected atoms to find the aromatic ring
        for neighbor in mol.GetAtomWithIdx(aryl_start).GetNeighbors():
            if neighbor.GetIsAromatic():
                # Check if the substituent forms an aromatic ring (at least 6-membered)
                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    if neighbor.GetIdx() in ring and len(ring) >= 6:
                        return True, "Contains 3-aryl-4H-chromen-4-one skeleton"
        return False, "Aryl substituent is not aromatic"

    # If no matches (unlikely due to earlier check), return False
    return False, "Unexpected error in structure matching"