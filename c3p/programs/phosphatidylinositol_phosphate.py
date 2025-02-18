"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: Phosphatidylinositol phosphate
A phosphatidylinositol phosphate (PIP, PIP2, PIP3, etc.) is a phosphoinositide that
contains a (myo-)inositol head group with one or more phosphate substituents on the ring,
a linking phosphate connecting the inositol to a glycerol backbone, and at least two acyl chains.
This version aims at reducing false positives (e.g. normal PI lipids) by requiring that the inositol head
has an extra phosphate substituent (other than the linking phosphate).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate (PIP/PIP2/PIP3) based on its SMILES string.
    
    The classifier checks that:
      0. The molecule does not contain any negative formal charges.
      1. It contains a myo-inositol head group (using a canonical inositol pattern).
      2. One of the inositol oxygens is connected to a 'linking phosphate' that in turn is attached (via an oxygen)
         to a carbon outside the inositol. This connects the inositol head to the glycerol backbone.
      3. At least one other inositol oxygen (i.e. not used for the linking phosphate) is phosphorylated. This extra
         phosphate on the ring differentiates phosphatidylinositol phosphates from the plain PI lipids.
      4. The molecule has at least two acyl ester groups (fatty acid chains attached via ester bonds) that are not
         part of any phosphate group.
      5. The molecular weight is above a minimal threshold typical for these lipids.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phosphatidylinositol phosphate, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 0. Reject molecules with any negative formal charges (to avoid deprotonated salts).
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() < 0:
            return False, "Molecule has negative formal charges; likely drawn as a salt form"
    
    # 1. Check for a myo-inositol head group using the canonical pattern.
    # This pattern corresponds to a fully hydroxylated cyclohexane: OC1C(O)C(O)C(O)C(O)C1O
    inositol = Chem.MolFromSmiles("OC1C(O)C(O)C(O)C(O)C1O")
    inositol_matches = mol.GetSubstructMatches(inositol)
    if not inositol_matches:
        return False, "Inositol head group not found"
    
    # Use the first matching inositol fragment.
    inositol_atom_indices = set(inositol_matches[0])
    
    # Variables to track the linking phosphate and extra phosphorylation on the inositol head.
    linking_found = False
    extra_phospho_found = False

    # We will loop over all atoms in the molecule whose index is in the inositol match.
    # For each inositol oxygen, examine if it is bound to phosphorus.
    for idx in inositol_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        # We're only interested in oxygen atoms from the inositol (the pattern gives both O and C,
        # but phosphate substitutions come via O atoms)
        if atom.GetSymbol() != "O":
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() != "P":
                continue
            # For each phosphorus neighbor, differentiate linking vs extra phosphorylation.
            # A linking phosphate is defined as: from the phosphorus (P_link)
            # there is an oxygen (other than the current one) that is bound to a carbon not in the inositol.
            is_linking = False
            for p_nbr in nbr.GetNeighbors():
                if p_nbr.GetIdx() == atom.GetIdx():
                    continue  # skip the inositol oxygen in question
                if p_nbr.GetSymbol() == "O":
                    # Check if this oxygen is attached to a carbon outside the inositol.
                    for second_nbr in p_nbr.GetNeighbors():
                        if second_nbr.GetSymbol() == "C" and second_nbr.GetIdx() not in inositol_atom_indices:
                            is_linking = True
                            break
                if is_linking:
                    break
            if is_linking:
                linking_found = True
            else:
                # If the phosphorus does not fulfill the linking role, assume it is a phosphoryl substituent on the inositol.
                extra_phospho_found = True

    if not linking_found:
        return False, "Linking phosphate connecting inositol to glycerol backbone not found"
    if not extra_phospho_found:
        return False, "No extra phosphate substituent on the inositol head found; likely a plain PI rather than PIP/PIP2/PIP3"
    
    # 4. Count acyl ester groups.
    # Look for ester linkages defined as an oxygen single-bonded to a carbon that is double-bonded to another oxygen.
    # Exclude any ester oxygen that is bonded to phosphorus (they might be part of a phosphate).
    acyl_count = 0
    seen_ester_bonds = set()
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify candidate pair: one atom must be oxygen and the other must be carbon.
        if a1.GetSymbol() == "O" and a2.GetSymbol() == "C":
            oxygen = a1
            carbon = a2
        elif a2.GetSymbol() == "O" and a1.GetSymbol() == "C":
            oxygen = a2
            carbon = a1
        else:
            continue
        
        # Exclude if this oxygen is attached to any phosphorus (it could be part of a phosphate group)
        if any(nb.GetSymbol() == "P" for nb in oxygen.GetNeighbors()):
            continue
        
        # Confirm the carbon is a carbonyl: it must have a double bond to an oxygen (other than our ester oxygen).
        carbonyl_found = False
        for nbr in carbon.GetNeighbors():
            # Skip the bond already considered.
            if nbr.GetIdx() == oxygen.GetIdx():
                continue
            if nbr.GetSymbol() != "O":
                continue
            bond_co = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
            if bond_co is not None and bond_co.GetBondType() == Chem.BondType.DOUBLE:
                carbonyl_found = True
                break
        if carbonyl_found:
            bond_id = tuple(sorted([oxygen.GetIdx(), carbon.GetIdx()]))
            if bond_id not in seen_ester_bonds:
                seen_ester_bonds.add(bond_id)
                acyl_count += 1

    if acyl_count < 2:
        return False, f"Found only {acyl_count} acyl ester group(s), need at least 2"
    
    # 5. Check molecular weight (typically these lipids are >500 Da).
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight too low ({mw:.1f} Da) for a phosphatidylinositol phosphate"
    
    return True, "Molecule contains a phosphorylated inositol head (with linking and extra phosphate), and at least two acyl chains"


# Example usage (for testing):
if __name__ == "__main__":
    # Example true positive (a PIP2 type molecule)
    test_smiles = "CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC"
    result, reason = is_phosphatidylinositol_phosphate(test_smiles)
    print(result, reason)