"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: Any organic heteromonocyclic compound based on a dihydropyrrole motif (pyrroline).
A pyrroline (here) is defined as a 5‐membered ring that contains exactly one nitrogen atom 
and exactly one double bond (non‐aromatic, internal to the ring).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_pyrroline(smiles: str):
    """
    Determines if the given molecule (via its SMILES string) contains a pyrroline (dihydropyrrole) core.
    The pyrroline motif is recognized as a 5-membered ring with exactly one nitrogen atom and exactly one double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a pyrroline core is found, else False.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings (as tuples of atom indices)
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    # Loop through each ring to test if it meets our criteria
    for ring in rings:
        if len(ring) != 5:
            continue  # Only consider 5-membered rings
        
        # Count the number of nitrogen atoms in the ring
        n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_nitrogen != 1:
            continue  # Not the desired heterocycle
        
        # Extract bonds which are strictly between the ring’s atoms
        ring_atom_set = set(ring)
        ring_bonds = []
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_atom_set and a2 in ring_atom_set:
                ring_bonds.append(bond)
                
        # For a simple 5-membered ring, we expect exactly 5 bonds.
        # Now count how many of these bonds are double bonds.
        # We check only for bonds with actual bond type "DOUBLE" (ignores aromatic bonds).
        double_bond_count = sum(1 for bond in ring_bonds if bond.GetBondType() == rdchem.BondType.DOUBLE)
        
        # We expect exactly one double bond in a dihydropyrrole (pyrroline) ring.
        if double_bond_count == 1:
            return True, "Found a 5-membered ring with one nitrogen and one double bond (dihydropyrrole core)."
    
    return False, "No 5-membered dihydropyrrole ring with one nitrogen and one double bond found."

# Example usage (this part can be omitted or used for testing):
if __name__ == '__main__':
    test_smiles = "CC(=O)C1=NCCC1"  # 2-acetyl-1-pyrroline example
    result, reason = is_pyrroline(test_smiles)
    print(result, reason)