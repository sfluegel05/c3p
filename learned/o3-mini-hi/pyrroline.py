"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: Any organic heteromonocyclic compound based on a dihydropyrrole motif (pyrroline).
A pyrroline (here) is defined as a 5‐membered ring (from the symmetric SSSR list) that contains
exactly one nitrogen atom and exactly one double bond (in a non‐aromatic, Kekulized ring).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_pyrroline(smiles: str):
    """
    Determines if the given molecule (via its SMILES string) contains a pyrroline (dihydropyrrole) core.
    The pyrroline motif is recognized as a 5-membered ring, when Kekulized, that has:
      - exactly one nitrogen atom in the ring, and
      - exactly one double bond among the ring bonds.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a qualifying pyrroline core is found, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Kekulize the molecule to convert aromatic bonds into explicit single/double bonds.
    # This helps make the double bond count in the ring explicit.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        # In some cases Kekulize may fail; we can continue with the original flags.
        pass

    # Get the set of small rings via the Symmetry-unique SSSR method:
    ssr = Chem.GetSymmSSSR(mol)
    if not ssr:
        return False, "No rings found in the molecule"
    
    # Loop over each 5-membered ring from the SSSR list
    for ring in ssr:
        if len(ring) != 5:
            continue  # we only consider rings with exactly 5 atoms
        
        # Count the number of nitrogen atoms in this ring.
        n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_nitrogen != 1:
            continue  # not the desired dihydropyrrole if there isn't exactly one nitrogen
        
        # Collect all bonds that connect two atoms in this ring.
        ring_atom_set = set(ring)
        ring_bonds = []
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_atom_set and a2 in ring_atom_set:
                ring_bonds.append(bond)
                
        # For a clean 5-membered ring we expect exactly 5 bonds.
        if len(ring_bonds) != 5:
            continue
        
        # Count the number of bonds that are explicit double bonds.
        double_bond_count = sum(1 for bond in ring_bonds if bond.GetBondType() == rdchem.BondType.DOUBLE)
        
        # In a dihydropyrrole (pyrroline) we expect exactly one double bond.
        if double_bond_count == 1:
            return True, "Found a 5-membered ring with one nitrogen and one double bond (dihydropyrrole core)."
    
    return False, "No 5-membered dihydropyrrole ring with one nitrogen and one double bond found."

# Example usage (for testing):
if __name__ == '__main__':
    # Test using one of the provided SMILES strings (2-acetyl-1-pyrroline)
    test_smiles = "CC(=O)C1=NCCC1"
    result, reason = is_pyrroline(test_smiles)
    print(result, reason)