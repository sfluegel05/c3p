"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
#!/usr/bin/env python3
"""
Classifies: 11β-hydroxy steroid
Defined as any 11-hydroxy steroid in which the hydroxy group at position 11 has beta configuration.
This heuristic routine first checks that the molecule has a steroid-like fused ring system (≥4 rings)
and then searches for a candidate beta-oriented hydroxyl group attached to a carbon in a 5 or 6-membered ring.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11β-hydroxy steroid based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES. If invalid, returns False.
      2. Checks if the molecule has a steroid-like fused ring system
         (heuristically at least 4 rings, many of which are 5- or 6-membered).
      3. Searches for a hydroxyl group – an oxygen atom attached to a carbon – 
         where the carbon is chiral and explicitly specified with '@@' (often used 
         to depict beta orientation in steroid SMILES) and is part of a 5- or 6-membered ring.
         
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as an 11β-hydroxy steroid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Heuristic: A steroid should have at least 4 fused rings.
    if len(rings) < 4:
        return False, f"Found only {len(rings)} rings; a steroid nucleus should have at least 4 fused rings"
        
    # To focus on the steroid core, we check if many of the rings are of size 5 or 6.
    ring_sizes = [len(r) for r in rings]
    num_56 = sum(1 for size in ring_sizes if size in (5,6))
    if num_56 < 3:
        return False, "Fewer than 3 rings of size 5 or 6; does not appear to have a typical steroid nucleus"
    
    # Search for a beta-oriented hydroxyl group.
    # We assume that the hydroxyl group should be a -OH attached to a chiral carbon,
    # whose stereochemistry is conveyed using the '@@' symbol (common way to represent beta in steroids).
    beta_hydroxy_found = False
    for atom in mol.GetAtoms():
        # Look for oxygen atoms (atomic number 8)
        if atom.GetAtomicNum() != 8:
            continue
        # Check that the oxygen is part of a hydroxyl group (i.e., connected to exactly one heavy atom)
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 1:
            continue
        carbon = neighbors[0]
        # Ensure the oxygen is a single-bonded -OH (and not e.g. part of a carbonyl)
        if mol.GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx()).GetBondTypeAsDouble() != 1:
            continue
        
        # Check if the attached carbon is chiral and carries a defined stereochemistry.
        # We then inspect the chiral tag: in RDKit, the stereochemistry is encoded as PARTY atom properties.
        if not carbon.HasProp('_CIPCode'):
            # If no CIP code is defined, we check the chiral tag manually.
            # The default chiral type is stored in the atom's chiralTag.
            if carbon.GetChiralTag() not in (Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
                                             Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW):
                continue
        
        # Retrieve the atom's chiral specification from its SMILES string.
        # One way is to re-generate the atom’s SMILES fragment.
        atom_smiles = Chem.MolToSmiles(Chem.MolFragmentToMol(mol, [carbon.GetIdx()]),
                                       isomericSmiles=True)
        # Look for the '@@' pattern as a heuristic for beta orientation.
        if '@@' not in atom_smiles:
            continue
        
        # Further, check that the carbon is part of one of the fused rings (of size 5 or 6)
        in_small_ring = False
        for ring in rings:
            if carbon.GetIdx() in ring and len(ring) in (5,6):
                in_small_ring = True
                break
        if not in_small_ring:
            continue
        
        # If we found such a candidate, assume it is the beta hydroxy group.
        beta_hydroxy_found = True
        break
        
    if not beta_hydroxy_found:
        return False, "No beta-oriented hydroxyl group found on a chiral carbon within a steroid ring"
    
    return True, "Steroid nucleus with at least 4 fused rings and a beta-oriented (11β) hydroxyl group detected"


# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: 11beta-hydroxyandrost-4-ene-3,17-dione
    test_smiles = "[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@]1(C)C(=O)CC[C@@]21[H]"
    result, reason = is_11beta_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)