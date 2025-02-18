"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: CHEBI:72025 flavone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    Flavones have a 2-aryl-1-benzopyran-4-one skeleton and substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the 2-arylchromen-4-one skeleton
    # Core is benzopyran-4-one with aryl group at position 2
    flavone_core = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)c(c2)-[a]")
    matches = mol.GetSubstructMatches(flavone_core)
    if not matches:
        return False, "No 2-arylchromen-4-one skeleton found"

    # Get ring information
    ring_info = mol.GetRingInfo()
    core_rings = set()
    for ring in ring_info.AtomRings():
        # Check if this ring is part of the core (benzopyran-4-one)
        # Core atoms are from the SMARTS match (indices 0-8 for the core + aryl connector)
        core_atoms = set(matches[0][:9])  # Adjust indices based on SMARTS
        if any(atom in core_atoms for atom in ring):
            core_rings.add(tuple(ring))

    # Check aryl group is a separate aromatic ring (not part of core)
    for match in matches:
        aryl_connector = match[-1]  # Last atom in match is the aryl connector
        for ring in ring_info.AtomRings():
            if aryl_connector in ring and tuple(ring) not in core_rings:
                # Check if the ring is aromatic and at least 6-membered
                is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
                if is_aromatic and len(ring) >= 6:
                    return True, "Contains 2-arylchromen-4-one skeleton with aromatic substituent"

    return False, "Aryl group not part of separate aromatic ring"