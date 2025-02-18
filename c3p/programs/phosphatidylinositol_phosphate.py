"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: Phosphatidylinositol phosphate
A phosphatidylinositol phosphate (PIP, PIP2, PIP3, etc.) is a phosphoinositide that,
in natural occurrence, is found as one of only seven compounds. They have a myo-inositol head
with a phosphate that links it to a glycerol backbone that carries two acyl (ester) chains,
plus additional phosphate(s) on the inositol ring.
This classifier uses a series of structural heuristics:
  0. Reject molecules with any negative formal charges.
  1. Ensure the molecule contains a myo-inositol head group.
  2. Count phosphorus atoms attached to the hydroxyls of this inositol. In a natural PIP family member,
     one of these is the “linking” phosphate (to the glycerol backbone) and there is at least one extra.
     We require that the total number of P attached to inositol (linker + extra) is 2, 3 or 4.
  3. Confirm that at least one of the P atoms shows linking behavior – with an oxygen that connects to a carbon
     not belonging to the inositol.
  4. Look for a glycerol backbone that is typical of diacyl lipids: two acyl ester bonds coming off the same carbon.
  5. Ensure the molecular weight is high enough typical for such lipids.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate (PIP/PIP2/PIP3)
    based on its SMILES string.
    
    The classifier applies the following tests:
      0. The molecule must not contain any negative formal charges.
      1. It must contain a myo-inositol head group.
      2. Count the phosphorus atoms attached to inositol's hydroxyl oxygens.
         (One is the linking phosphate to the glycerol and the others are extra phosphates.)
         The total must be 2, 3, or 4.
      3. At least one phosphorus (the linker) must have a neighbor oxygen that in turn is attached
         to a carbon not part of the inositol.
      4. The molecule must contain a glycerol backbone with two acyl ester (fatty acid) chains. In our
         heuristic we require that the two acyl ester bonds originate from the same carbon (the typical
         diacyl glycerol).
      5. The molecular weight must be high (here, above 500 Da).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phosphatidylinositol phosphate, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 0. Reject molecules with negative formal charges (to skip salt forms)
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() < 0:
            return False, "Molecule has negative formal charges; likely drawn as a salt form"

    # 1. Identify a myo-inositol head.
    # We use a canonical myo-inositol: note the stereochemistry is set to match typical myo-inositol.
    inositol = Chem.MolFromSmiles("O[C@@H]1[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O")
    inositol_matches = mol.GetSubstructMatches(inositol)
    if not inositol_matches:
        return False, "Inositol head group not found"
    # Use the first matching substructure for further analysis.
    inositol_atom_indices = set(inositol_matches[0])

    # 2. Count phosphorus atoms attached directly to the inositol oxygens.
    # (They will be attached to the –OH groups of the inositol).
    p_attached = {}  # phosphorus idx -> list of inositol oxygen indices attached to it
    for idx in inositol_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != "O":
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() == "P":
                p_idx = nbr.GetIdx()
                p_attached.setdefault(p_idx, []).append(idx)
    num_inositol_ps = len(p_attached)
    # We expect a linking phosphate plus at least one extra phosphate.
    # Thus, total P count on inositol should be 2, 3 or 4.
    if num_inositol_ps not in (2, 3, 4):
        return False, ("Found %d phosphate(s) attached to the inositol head; "
                       "expected 2, 3, or 4 (linking plus extra phosphate(s))." % num_inositol_ps)

    # 3. Check for a linking phosphate functionality.
    # For at least one phosphorus, one of its attached oxygens (that is not part of inositol)
    # should be bound to a carbon not in the inositol fragment.
    linking_found = False
    for p_idx, o_indices in p_attached.items():
        p_atom = mol.GetAtomWithIdx(p_idx)
        for o_neigh in p_atom.GetNeighbors():
            if o_neigh.GetIdx() in inositol_atom_indices:
                continue
            if o_neigh.GetSymbol() != "O":
                continue
            # Look for attachment to a non-inositol carbon.
            for second in o_neigh.GetNeighbors():
                if second.GetSymbol() == "C" and second.GetIdx() not in inositol_atom_indices:
                    linking_found = True
                    break
            if linking_found:
                break
        if linking_found:
            break
    if not linking_found:
        return False, "Linking phosphate (connecting inositol head to glycerol backbone) not found"

    # 4. Count acyl ester bonds and check for a glycerol backbone.
    # We define an acyl ester bond as a single bond between an oxygen and a carbon where:
    #   - The oxygen is not attached to any phosphorus.
    #   - The carbon (on the acyl side) is double-bonded to an oxygen (carbonyl).
    #
    # To ensure these acyl groups come from a common glycerol backbone, we record, for each carbon,
    # how many such ester bonds it forms.
    acyl_bonds_by_carbon = {}
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify a candidate pair: one atom must be oxygen, the other carbon.
        if a1.GetSymbol() == "O" and a2.GetSymbol() == "C":
            oxygen = a1
            carbon = a2
        elif a2.GetSymbol() == "O" and a1.GetSymbol() == "C":
            oxygen = a2
            carbon = a1
        else:
            continue
        
        # Exclude oxygen if it is attached to any phosphorus (likely part of a phosphate)
        if any(nb.GetSymbol() == "P" for nb in oxygen.GetNeighbors()):
            continue
        
        # Confirm that the carbon is part of a carbonyl:
        carbonyl_found = False
        for nbr in carbon.GetNeighbors():
            # Skip the oxygen in the ester bond
            if nbr.GetIdx() == oxygen.GetIdx():
                continue
            if nbr.GetSymbol() != "O":
                continue
            bond_co = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
            if bond_co is not None and bond_co.GetBondType() == Chem.BondType.DOUBLE:
                carbonyl_found = True
                break
        if carbonyl_found:
            cid = carbon.GetIdx()
            acyl_bonds_by_carbon[cid] = acyl_bonds_by_carbon.get(cid, 0) + 1

    # Check if any carbon has two or more such ester bonds (indicative of a glycerol with two acyl chains).
    if not acyl_bonds_by_carbon or max(acyl_bonds_by_carbon.values()) < 2:
        return False, "No glycerol backbone with at least two acyl ester groups detected (acyl chains not on the same carbon)"

    # 5. Check molecular weight (a typical phosphatidylinositol phosphate is heavy).
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, "Molecular weight too low (%.1f Da) for a phosphatidylinositol phosphate" % mw

    return True, ("Molecule contains a myo-inositol head with the expected phosphate substitutions "
                  "and a glycerol backbone with two acyl chains")

# Example usage (for testing):
if __name__ == "__main__":
    # Test with a known true positive:
    test_smiles = "CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCCC"
    result, reason = is_phosphatidylinositol_phosphate(test_smiles)
    print(result, reason)