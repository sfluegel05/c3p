"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Organochlorine compounds that are polychlorinated dibenzodioxines and related compounds.
Definition:
  - Must contain at least 2 halogen atoms (Cl or Br).
  - Molecular weight must be <=700 Da.
  - Must have one of the following aromatic scaffolds:
         • Biphenyl (two benzene rings connected by a single bond)
         • Diaryl ether (two benzene rings connected by an oxygen atom)
         • Dibenzodioxin
         • Dibenzofuran
These criteria aim to capture the core structure of polychlorinated persistent organic pollutants.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated 
    dibenzodioxines and related compounds.

    Criteria:
      1. Contains at least 2 halogen atoms (Cl (Z=17) or Br (Z=35)).
      2. Its molecular weight is <= 700 Da.
      3. Contains one or more characteristic aromatic scaffolds:
            - Biphenyl (two benzene rings joined by a single bond)
            - Diaryl ether (two benzene rings joined by an oxygen)
            - Dibenzodioxin: pattern c1ccc2Oc3ccccc3O2c1
            - Dibenzofuran: pattern c1ccc2oc3ccccc3c2c1

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the decision.
    """
    
    # Convert SMILES to RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: At least 2 halogen atoms (chlorine or bromine).
    halogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in (17, 35)]
    if len(halogen_atoms) < 2:
        return False, "Insufficient halogen atoms (need at least 2 Cl or Br substituents)"
    
    # Criterion 2: Molecular weight must be <= 700 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 700:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da); should be <= 700 Da"
    
    # Criterion 3: Check for aromatic scaffolds.
    # Define SMARTS for biphenyl (directly connected benzene rings)
    biphenyl_smarts = "c1ccccc1-c2ccccc2"
    # Diaryl ether: two benzene rings connected through an oxygen. This is often present in aliased structures.
    diaryl_ether_smarts = "c1ccccc1Oc2ccccc2"
    # Dibenzodioxin.
    dibenzodioxin_smarts = "c1ccc2Oc3ccccc3O2c1"
    # Dibenzofuran.
    dibenzofuran_smarts = "c1ccc2oc3ccccc3c2c1"
    
    # Compile SMARTS patterns.
    biphenyl_pattern = Chem.MolFromSmarts(biphenyl_smarts)
    diaryl_ether_pattern = Chem.MolFromSmarts(diaryl_ether_smarts)
    dibenzodioxin_pattern = Chem.MolFromSmarts(dibenzodioxin_smarts)
    dibenzofuran_pattern = Chem.MolFromSmarts(dibenzofuran_smarts)
    
    scaffolds_found = []
    if mol.HasSubstructMatch(biphenyl_pattern):
        scaffolds_found.append("biphenyl")
    if mol.HasSubstructMatch(diaryl_ether_pattern):
        scaffolds_found.append("diaryl ether")
    if mol.HasSubstructMatch(dibenzodioxin_pattern):
        scaffolds_found.append("dibenzodioxin")
    if mol.HasSubstructMatch(dibenzofuran_pattern):
        scaffolds_found.append("dibenzofuran")
    
    if not scaffolds_found:
        return False, "No recognized aromatic scaffold (biphenyl, diaryl ether, dibenzodioxin, or dibenzofuran) found"
    
    reason = (
        f"Molecule has a {' and '.join(scaffolds_found)} scaffold, "
        f"{len(halogen_atoms)} halogen substituents, and a molecular weight of {mol_wt:.1f} Da"
    )
    
    return True, reason

# Example usage (for testing):
if __name__ == '__main__':
    # Test with 4,4'-dichlorobiphenyl (should be True).
    test_smiles = "Clc1ccc(cc1)-c1ccc(Cl)cc1"
    result, msg = is_polychlorinated_dibenzodioxines_and_related_compounds(test_smiles)
    print(result, msg)