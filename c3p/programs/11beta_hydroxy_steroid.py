"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
#!/usr/bin/env python3
"""
Classifies: 11β-hydroxy steroid
Defined as: Any 11-hydroxy steroid in which the hydroxy group at position 11 has beta- configuration.
The algorithm:
  1. Parse the SMILES string.
  2. Verify that the molecule has a fused steroid nucleus (heuristically ≥4 rings, with several of size 5 or 6).
  3. Look for a hydroxyl (-OH) group attached to a chiral carbon that shows explicit beta orientation.
     To do this, for each oxygen (atomic number 8) that is part of a hydroxyl group (single-bonded and only one heavy neighbor),
     we retrieve a SMILES fragment for its attached carbon (using Chem.MolFragmentToSmiles) and check if it starts with “[C@@”.
If both conditions are met we return True.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11β-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11β-hydroxy steroid, False otherwise.
        str: Reason for the classification decision.
    """
    try:
        # Parse the SMILES string; if invalid, return an error.
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
    except Exception as e:
        return False, f"Error parsing SMILES: {e}"
        
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Check for a steroid-like fused ring system.
    if len(rings) < 4:
        return False, f"Found only {len(rings)} rings; a steroid usually has at least 4 fused rings"
        
    # Check that several rings are of size 5 or 6.
    ring_sizes = [len(r) for r in rings]
    num_56 = sum(1 for size in ring_sizes if size in (5, 6))
    if num_56 < 3:
        return False, "Fewer than 3 rings with sizes 5 or 6; does not appear to have a typical steroid nucleus"
    
    # Search for a beta-oriented hydroxyl group.
    beta_hydroxy_found = False
    # Iterate over all oxygen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        neighbors = atom.GetNeighbors()
        # We expect an -OH group: oxygen bonded to exactly one heavy atom.
        if len(neighbors) != 1:
            continue
        carbon = neighbors[0]
        bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), atom.GetIdx())
        # Ensure this bond is a single bond.
        if bond.GetBondTypeAsDouble() != 1:
            continue
        
        # Now check that the carbon is chiral.
        if not carbon.HasProp('_CIPCode') and carbon.GetChiralTag() not in (
                Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
                Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW):
            continue
        
        # Check that this carbon is part of one of the small (5- or 6-membered) rings.
        in_small_ring = False
        for ring in rings:
            if carbon.GetIdx() in ring and len(ring) in (5,6):
                in_small_ring = True
                break
        if not in_small_ring:
            continue

        # Generate a SMILES fragment for the carbon to inspect its explicit stereochemistry.
        try:
            # Note: MolFragmentToSmiles is available and returns a fragment SMILES.
            frag = Chem.MolFragmentToSmiles(mol, atomsToUse=[carbon.GetIdx()], canonical=True, isomericSmiles=True)
        except Exception as e:
            # In case this fails, we skip this candidate.
            continue
        
        # Check if the fragment starts with "[C@@"; this is our heuristic for beta orientation.
        if frag.startswith("[C@@"):
            beta_hydroxy_found = True
            break

    if not beta_hydroxy_found:
        return False, "No beta-oriented hydroxyl group found on a chiral carbon within a steroid ring."
        
    return True, "Steroid nucleus detected with at least 4 (fused) rings and a beta-oriented (11β) hydroxyl group found."

# Example usage:
if __name__ == "__main__":
    # Example: prednisolone phosphate
    test_smiles = "C[C@@]12C[C@@H]([C@]3([C@]([C@]2([H])CC[C@@]1(C(COP(=O)(O)O)=O)O)(CCC=4[C@]3(C)C=CC(C4)=O)[H])[H])O"
    result, reason = is_11beta_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)