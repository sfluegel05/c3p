"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: Sesquiterpenoid
Definition:
  Any terpenoid derived from a sesquiterpene (C15 skeleton) even when rearranged or modified.
  Here we require that the overall molecule:
    - Contains only C, H, and O atoms (no aromatic rings);
    - Has 12–18 carbon atoms (approximately the 15 carbons expected);
    - Has a low molecular weight (150–400 Da);
    - Has a degree of unsaturation (DBE, calculated as C - H/2 + 1) between 3 and 7.
Note: This heuristic is not perfect, but aims to improve on previous scaffold‐based counts.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are derived from sesquiterpenes (starting from a C15 skeleton)
    even when rearranged. The classification is done by applying several filters:
      (1) only C, H, O are allowed;
      (2) no aromatic atoms (sesquiterpenoid cores are usually aliphatic);
      (3) total carbon count between 12 and 18;
      (4) molecular weight between 150 and 400 Da;
      (5) degree of unsaturation (DBE = C - H/2 + 1) is between 3 and 7.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if classified as a sesquiterpenoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check that the molecule contains only allowed elements: C (6), H (1), and O (8)
    allowed_atomic_numbers = {1, 6, 8}  # only H, C, and O allowed
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_numbers:
            return False, f"Contains atom {atom.GetSymbol()}, not allowed for a typical sesquiterpenoid"
    
    # 2. Check that the molecule does not have aromatic atoms.
    # (Terpenoid skeletons are usually aliphatic; if any atom is aromatic, likely not a sesquiterpenoid.)
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if aromatic_atoms:
        return False, "Contains aromatic atoms, which is not typical for a sesquiterpenoid structure"
    
    # 3. Count carbons and hydrogens in the entire molecule.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if not (12 <= c_count <= 18):
        return False, f"Molecule contains {c_count} carbons, not consistent with expected ~15 (range 12–18) for a sesquiterpene-derived structure"
    
    # 4. Check molecular weight (should be in a low range typical for sesquiterpenoids).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (150 <= mol_wt <= 400):
        return False, f"Molecular weight {mol_wt:.1f} Da not in expected range (150–400 Da) for sesquiterpenoids"
    
    # 5. Estimate degree of unsaturation (DBE) using the formula:
    #    DBE = C - H/2 + 1. (Since only C, H, and O are present, and O does not affect DBE)
    dbe = c_count - (h_count / 2) + 1
    # Allow a tolerance in DBE (terpenoids typically have 3–7 degrees of unsaturation)
    if not (3 <= dbe <= 7):
        return False, f"Degree of unsaturation (DBE) calculated as {dbe:.1f} is not in the expected range (3–7) for sesquiterpenoids"
    
    return True, f"Molecule has {c_count} carbons, MW {mol_wt:.1f} Da, and DBE {dbe:.1f}, consistent with a sesquiterpene-derived skeleton"

# Example usage (for testing purposes)
if __name__ == '__main__':
    # Use one of the provided examples, e.g. indicanone:
    test_smiles = "CC(=C)[C@@H]1CCC(CO)=C2CC(=O)C(C)=C2C1"
    result, reason = is_sesquiterpenoid(test_smiles)
    print(result, reason)