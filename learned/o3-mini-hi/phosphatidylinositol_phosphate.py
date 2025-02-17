"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: Phosphatidylinositol phosphate 
A phosphatidylinositol phosphate (PIP, PIP2, etc.) is expected to contain:
  - a (myo-)inositol head group,
  - a linking phosphate connecting the inositol to a glycerol backbone,
  - at least two acyl chains attached via ester bonds to the glycerol.
In this version we also reject molecules having any negative formal charge
to avoid misclassifying salt forms.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    The classifier checks that:
      1. The molecule does not contain any negative formal charges.
      2. It contains an inositol head group.
      3. At least one hydroxyl on the inositol is linked (through an oxygen) to a phosphate 
         that in turn is attached to a carbon (as expected for the linking phosphate).
      4. It has at least two acyl ester groups (fatty acid chains attached via ester bonds)
         that are not part of a phosphate.
      5. Its molecular weight is above a minimal threshold typical for these lipids.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a phosphatidylinositol phosphate, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 0. Reject molecules with any negative formal charges (to remove deprotonated(salt) forms).
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() < 0:
            return False, "Molecule has negative formal charges; likely drawn as a salt form"
    
    # 1. Check for an inositol head group.
    # The canonical myo-inositol (with hydroxyls) is represented by:
    inositol = Chem.MolFromSmiles("OC1C(O)C(O)C(O)C(O)C1O")
    inositol_matches = mol.GetSubstructMatches(inositol)
    if not inositol_matches:
        return False, "Inositol head group not found"
    
    # 2. Look for a linking phosphate that connects the inositol head group to a glycerol-like carbon.
    # The idea is: one of the hydroxyl oxygens of the inositol should be bonded to a phosphorus
    # that in turn is bonded (through another oxygen) to a carbon (from the glycerol backbone).
    linking_found = False
    inositol_atom_indices = set(inositol_matches[0])
    # Iterate over the atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        # If this oxygen is (likely) from the inositol head, it should be in the matching set.
        if atom.GetIdx() not in inositol_atom_indices:
            continue
        # Check if this inositol oxygen is bonded to phosphorus.
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() == "P":
                # Now from the phosphorus, check if one of its other neighbors is an oxygen
                # that in turn is bound to a carbon.
                for p_nbr in nbr.GetNeighbors():
                    if p_nbr.GetIdx() == atom.GetIdx():
                        continue
                    if p_nbr.GetSymbol() == "O":
                        for second_nbr in p_nbr.GetNeighbors():
                            if second_nbr.GetSymbol() == "C":
                                linking_found = True
                                break
                    if linking_found:
                        break
                if linking_found:
                    break
        if linking_found:
            break
    if not linking_found:
        return False, "Linking phosphate connecting inositol to glycerol backbone not found"
    
    # 3. Count acyl ester groups.
    # We search for an ester linkage defined as: an oxygen single-bonded to a carbon which is double-bonded to an oxygen.
    # And we skip any oxygen that is attached to a phosphorus (to avoid counting the linking phosphate).
    acyl_count = 0
    seen_esters = set()
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify pair where one is oxygen (ester oxygen candidate) and the other is carbon.
        if a1.GetSymbol() == "O" and a2.GetSymbol() == "C":
            oxygen = a1
            carbon = a2
        elif a2.GetSymbol() == "O" and a1.GetSymbol() == "C":
            oxygen = a2
            carbon = a1
        else:
            continue
        
        # Exclude if this oxygen is bonded to any phosphorus (could be part of the phosphate linking group).
        if any(nb.GetSymbol() == "P" for nb in oxygen.GetNeighbors()):
            continue
        
        # Confirm that the carbon is a carbonyl carbon (has a double bond to an oxygen).
        has_carbonyl = False
        for nbr in carbon.GetNeighbors():
            if nbr.GetSymbol() != "O" or nbr.GetIdx() == oxygen.GetIdx():
                continue
            bond_co = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
            if bond_co is not None and bond_co.GetBondType() == Chem.BondType.DOUBLE:
                has_carbonyl = True
                break
        if has_carbonyl:
            bond_id = tuple(sorted([oxygen.GetIdx(), carbon.GetIdx()]))
            if bond_id not in seen_esters:
                seen_esters.add(bond_id)
                acyl_count += 1
    
    if acyl_count < 2:
        return False, f"Found only {acyl_count} acyl ester group(s), need at least 2"
    
    # 4. Check molecular weight (typically these lipids are >500 Da).
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight too low ({mw:.1f} Da) for a phosphatidylinositol phosphate"
    
    return True, "Molecule contains a non-charged inositol head, a proper linking phosphate, and at least two acyl chains"

# Example usage (for testing):
if __name__ == "__main__":
    # Test on one true positive (should return True)
    test_smiles = "CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC"
    result, reason = is_phosphatidylinositol_phosphate(test_smiles)
    print(result, reason)