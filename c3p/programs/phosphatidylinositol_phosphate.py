"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: Phosphatidylinositol phosphate (a phosphoinositide, one of the seven naturally occurring).
A phosphatidylinositol phosphate (PIP, PIP2, etc.) is expected to contain:
  - a (myo-)inositol head group,
  - a phosphate ester linking the inositol head to a glycerol backbone (and inositol itself may be additionally phosphorylated),
  - at least two acyl chains attached via ester bonds to the glycerol.
This classifier uses a combination of substructure matching and bond-level analysis.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    The classifier checks for the presence of a myo-inositol head group, a linking phosphate (with additional phosphorylation on inositol),
    and at least two acyl ester (fatty acid) chains.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a phosphatidylinositol phosphate, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for an inositol head group.
    # We use the drawn myo-inositol SMILES (a cyclohexane with hydroxyl substituents).
    inositol = Chem.MolFromSmiles("OC1C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(inositol):
        return False, "Inositol head group not found"
    
    # 2. Ensure the molecule is a phosphorylated species (i.e. a phosphate ester is present linking inositol).
    # For phosphatidylinositol phosphates, we expect at least one additional phosphate on the head group,
    # so the total phosphorus count should be at least 2.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) < 2:
        return False, "Not enough phosphorus atoms; likely not a phosphorylated inositol (phosphatidylinositol phosphate)"
    
    # 3. Count acyl ester groups (i.e. the fatty acid chains attached to the glycerol backbone)
    # We'll count an ester group by scanning bonds: if an oxygen is bonded to a carbon that has a double-bond
    # to another oxygen (a carbonyl) then we consider that an acyl ester bond.
    # We also ensure that the oxygen is not part of a phosphate (i.e. not bonded to any phosphorus).
    acyl_count = 0
    seen_esters = set()  # To avoid double-counting bonds
    for bond in mol.GetBonds():
        # Only consider single bonds (the ester linkage is a single bond: R-O-C(=O)R')
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        
        # Check if one atom is oxygen and the other is carbon.
        if a1.GetSymbol() == 'O' and a2.GetSymbol() == 'C':
            oxygen = a1
            carbon = a2
        elif a2.GetSymbol() == 'O' and a1.GetSymbol() == 'C':
            oxygen = a2
            carbon = a1
        else:
            continue
        
        # Exclude oxygens that are part of a phosphate group (i.e. bonded to phosphorus).
        if any(nb.GetSymbol() == 'P' for nb in oxygen.GetNeighbors()):
            continue
        
        # Check if the carbon is a carbonyl carbon:
        # It must have at least one double bond to an oxygen (other than the oxygen we just saw).
        has_carbonyl = False
        for nb in carbon.GetNeighbors():
            if nb.GetSymbol() != 'O' or nb.GetIdx() == oxygen.GetIdx():
                continue
            bond_co = mol.GetBondBetweenAtoms(carbon.GetIdx(), nb.GetIdx())
            if bond_co and bond_co.GetBondType() == Chem.BondType.DOUBLE:
                has_carbonyl = True
                break
        
        if has_carbonyl:
            # Count this ester linkage if not already counted.
            bond_id = tuple(sorted([oxygen.GetIdx(), carbon.GetIdx()]))
            if bond_id not in seen_esters:
                seen_esters.add(bond_id)
                acyl_count += 1
                
    if acyl_count < 2:
        return False, f"Found only {acyl_count} acyl ester group(s), need at least 2"
    
    # 4. Optional extra check: the overall molecular weight should be in the range expected 
    # for these lipids (typically >500 Da).
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight too low ({mw:.1f} Da) for a phosphatidylinositol phosphate"
    
    return True, "Molecule contains a phosphorylated inositol head, a linking phosphate, and at least two acyl chains"

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = "CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H]1O)OC(=O)CCCCCCC"
    result, reason = is_phosphatidylinositol_phosphate(test_smiles)
    print(result, reason)