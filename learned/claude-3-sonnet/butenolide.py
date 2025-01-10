"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: CHEBI:35631 butenolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone containing a 2-furanone skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Multiple SMARTS patterns to catch different butenolide variations
    patterns = [
        # Basic 2-furanone pattern with any substitution
        "O=C1OCC=C1*",
        # Alternative pattern with double bond in different position
        "O=C1OC=CC1*",
        # More general pattern for substituted butenolides
        "O=C1OC([*,H])=C([*,H])C1([*,H])",
        # Pattern for highly substituted variants
        "O=C1OC([*,H])([*,H])C([*,H])=C1([*,H])",
        # Pattern catching some complex cases
        "O=C1OC([*,H])C=C1",
        # Pattern for fused ring systems
        "O=C1OCC2=C1[*]2"
    ]

    # Check each pattern
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            # Verify it's part of a 5-membered ring
            matches = mol.GetSubstructMatches(patt)
            ring_info = mol.GetRingInfo()
            for match in matches:
                # Check if matched atoms are part of a 5-membered ring
                for ring in ring_info.AtomRings():
                    if any(idx in ring for idx in match):
                        if len(ring) == 5:
                            # Additional verification: check for lactone group
                            lactone_atoms = set(match)
                            ring_atoms = set(ring)
                            # Verify oxygen and carbonyl carbon are in the ring
                            for atom_idx in lactone_atoms.intersection(ring_atoms):
                                atom = mol.GetAtomWithIdx(atom_idx)
                                if atom.GetAtomicNum() == 8:  # Oxygen
                                    neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
                                    if 6 in neighbors:  # Connected to carbon
                                        return True, "Contains butenolide ring system"

    # If no valid butenolide found
    return False, "No butenolide ring system found"