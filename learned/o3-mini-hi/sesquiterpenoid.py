"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: Sesquiterpenoid
Definition:
  Any terpenoid derived from a sesquiterpene (C15 skeleton) even if the skeleton
  has been rearranged or modified (e.g., by removal of one or more skeletal atoms).
  For our heuristic, we require:
    - Only C, H, and O atoms are allowed.
    - No aromatic atoms (sesquiterpene cores are usually aliphatic).
    - Total carbon count roughly between 12 and 21.
    - Molecular weight between 150 and 500 Da.
    - As an extra check, if the molecule is completely acyclic and bears a carboxylic acid,
      it is likely a fatty acid rather than a sesquiterpenoid.
Note: This heuristic is approximate and may not cover every edge case.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    
    The molecule must:
      (1) Contain only C, H, and O atoms.
      (2) Not contain any aromatic atoms.
      (3) Have a carbon count between 12 and 21.
      (4) Have a molecular weight between 150 and 500 Da.
      (5) (Heuristic) If acyclic and contains a carboxylic acid group, then it is likely a fatty acid.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a sesquiterpenoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check allowed elements: only hydrogen, carbon, and oxygen.
    allowed_atomic_nums = {1, 6, 8}  # H, C, and O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()}, not allowed for a typical sesquiterpenoid"
    
    # 2. Reject if any atom is aromatic.
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if aromatic_atoms:
        return False, "Contains aromatic atoms, which is not typical for a sesquiterpenoid structure"
    
    # 3. Count carbon atoms.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (12 <= c_count <= 21):
        return False, f"Molecule contains {c_count} carbons, not consistent with a typical sesquiterpene-derived structure (expected roughly 12–21 carbons)"
    
    # 4. Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (150 <= mol_wt <= 500):
        return False, f"Molecular weight {mol_wt:.1f} Da not in expected range (150–500 Da) for sesquiterpenoids"
    
    # 5. (Heuristic) If the molecule is acyclic and contains a carboxylic acid group,
    #    then it is likely a fatty acid.
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        # Check for a carboxylic acid functional group.
        # SMARTS pattern for a carboxylic acid group: carbonyl (C(=O)) bonded to -OH.
        acid_smarts = Chem.MolFromSmarts("[CX3](=O)[OH]")
        if mol.HasSubstructMatch(acid_smarts):
            return False, "Acyclic structure with a carboxylic acid group – likely a fatty acid rather than a sesquiterpenoid"
    
    return True, f"Molecule has {c_count} carbons, MW {mol_wt:.1f} Da, and passes allowed element and aliphatic checks, consistent with a sesquiterpene-derived structure"
    
# Example usage (for testing purposes)
if __name__ == '__main__':
    # You can test with one of the provided examples, e.g. indicanone:
    test_smiles = "CC(=C)[C@@H]1CCC(CO)=C2CC(=O)C(C)=C2C1"
    result, reason = is_sesquiterpenoid(test_smiles)
    print(result, reason)