"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is a five-membered cyclic ester (lactone ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a five-membered ring containing an ester group (O connected to C in a ring with C=O)
    # The pattern matches any oxygen in a 5-membered ring where there's a carbonyl group in the same ring
    gamma_lactone_pattern = Chem.MolFromSmarts("[O;r5].[C;r5]=O")
    
    # Check for matches with the pattern
    matches = mol.GetSubstructMatches(gamma_lactone_pattern)
    
    if not matches:
        return False, "No five-membered lactone ring found"
    
    # Check each possible pair of O and C=O in the same ring
    ring_info = mol.GetRingInfo()
    for match in matches:
        o_atom = match[0]
        c_atom = match[1]
        # Find rings containing both atoms
        common_rings = [r for r in ring_info.AtomRings() if o_atom in r and c_atom in r]
        for ring in common_rings:
            if len(ring) == 5:
                # Check if the oxygen is bonded to a carbon in the ring
                o_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(o_atom).GetNeighbors()]
                # Ensure the oxygen is bonded to a carbon in the ring
                if any(n in ring for n in o_neighbors):
                    # Check if the carbonyl carbon is part of an ester (bonded to an oxygen outside the ring)
                    # The carbonyl oxygen should not be in the ring
                    carbonyl_o = None
                    for bond in mol.GetAtomWithIdx(c_atom).GetBonds():
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            other = bond.GetOtherAtomIdx(c_atom)
                            if mol.GetAtomWithIdx(other).GetAtomicNum() == 8 and other not in ring:
                                carbonyl_o = other
                                break
                    if carbonyl_o is not None:
                        return True, "Contains a five-membered lactone ring"
    
    return False, "No valid five-membered lactone ring found"