"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
Definition: Any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6.
A heuristic approach is used:
  1. Look for at least one six-membered ring that contains a C=C double bond (as a proxy for the Δ5 double bond).
  2. Look for a chiral carbon in a ring that bears a hydroxyl group (proxy for 3β-hydroxy).
Note: Due to the complex nature of steroid nomenclature, this method uses a simplified approach.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if the given SMILES string represents a 3β-hydroxy-Δ5-steroid.
    
    The method uses two heuristic checks:
      1. The molecule must contain at least one six-membered ring that has a double bond.
         This approximates the Δ5 (i.e. C5=C6) feature found in many steroids.
      2. The molecule must contain at least one ring-bound chiral carbon that has an attached hydroxyl group.
         This serves as a proxy for the 3β-hydroxy group.

    Args:
        smiles (str): SMILES string of a molecule.
        
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Brief reason for the classification.
    """
    # Parse the SMILES to create a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic 1: Look for a six-membered ring that has a double bond.
    ring_with_db_found = False
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:  # six-membered ring
            # Check each bond in the ring (remembering the ring is cyclic)
            for i in range(len(ring)):
                j = ring[(i + 1) % len(ring)]
                bond = mol.GetBondBetweenAtoms(ring[i], j)
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    ring_with_db_found = True
                    break
            if ring_with_db_found:
                break
    if not ring_with_db_found:
        return False, "No six-membered ring with a double bond (Δ5 feature) found"
    
    # Heuristic 2: Look for a chiral carbon in a ring that carries a hydroxyl group.
    # This is our proxy for a 3β-hydroxy substituent.
    hydroxy_found = False
    for atom in mol.GetAtoms():
        # Consider carbon atoms in any ring
        if atom.GetAtomicNum() == 6 and atom.IsInRing():
            # Check if stereochemistry is set (as a heuristic for beta orientation)
            # _CIPCode is usually set if the stereochemistry is assigned.
            if atom.HasProp('_CIPCode'):
                # Look through the neighbors for an oxygen (typically part of an -OH group)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8:
                        # Check that the oxygen is in a hydroxyl (i.e. has a hydrogen attached)
                        # GetTotalNumHs() tells the number of attached hydrogens.
                        if neighbor.GetTotalNumHs() >= 1:
                            hydroxy_found = True
                            break
        if hydroxy_found:
            break
    if not hydroxy_found:
        return False, "No ring-bound chiral carbon with an attached hydroxyl (3β-hydroxy) found"
    
    return True, ("Molecule contains at least one six-membered ring with a double bond (Δ5 feature) and "
                  "a ring-bound chiral carbon bearing a hydroxyl group (proxy for 3β-hydroxy), "
                  "which is consistent with a 3β-hydroxy-Δ5-steroid.")