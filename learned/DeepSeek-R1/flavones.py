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
    # Chromen-4-one core (benzopyran-4-one) with aryl group at position 2
    # Core: benzene fused to pyrone (O=C1C=CC2=C1C=CC=C2)
    # Aryl group attached to C2 of pyrone (atom index 4 in SMARTS)
    flavone_core = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)cc2-[a]")
    matches = mol.GetSubstructMatches(flavone_core)
    if not matches:
        return False, "No 2-arylchromen-4-one skeleton found"

    # Get ring information
    ring_info = mol.GetRingInfo()
    core_rings = set()
    for ring in ring_info.AtomRings():
        if any(atom_idx in matches[0][:6] for atom_idx in ring):  # Core atoms are first 6 in SMARTS pattern
            core_rings.add(ring)

    # Check aryl group is a separate aromatic ring
    for match in matches:
        aryl_connector = match[-1]  # Last atom in match is the aryl connector
        # Check if aryl connector is part of a separate aromatic ring
        for ring in ring_info.AtomRings():
            if aryl_connector in ring and ring not in core_rings:
                # Check if the ring is aromatic
                is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
                if is_aromatic and len(ring) >= 6:  # At least benzene-like
                    return True, "Contains 2-arylchromen-4-one skeleton with aromatic substituent"

    return False, "Aryl group not part of separate aromatic ring"