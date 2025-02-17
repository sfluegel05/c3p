"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Organochlorine compounds that are polychlorinated dibenzodioxines and related compounds.
These include polychlorinated dibenzofurans as well as polychlorinated (or polybrominated) biphenyls.
The criteria in this revised version are:
  - The molecule must have at least 2 halogen atoms (Cl or Br).
  - Its molecular weight is capped at 700 Da (expanded from the original 600 Da cutoff).
  - It must contain at least one of the characteristic aromatic scaffolds:
         (a) biphenyl,
         (b) dibenzodioxin,
         (c) dibenzofuran.
  Note: We no longer require exactly 2 aromatic rings so that extra substituents do not lead to misclassification.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxines and related compounds.
    The classification applies the following criteria:
      1. The molecule must contain at least 2 halogen atoms (chlorine or bromine).
      2. The molecular weight is required to be below or equal to 700 Da (to allow slightly larger, but still persistent, pollutants).
      3. The molecule must contain one of the characteristic aromatic scaffolds:
             - biphenyl,
             - dibenzodioxin,
             - dibenzofuran.
         We use recursive SMARTS patterns to find these scaffolds regardless of substituents.
        
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count halogen atoms: chlorine (atomic number 17) and bromine (atomic number 35)
    halogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in (17, 35)]
    if len(halogen_atoms) < 2:
        return False, "Insufficient halogen atoms (need at least 2 Cl or Br substituents)"
    
    # Compute molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 700:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da); typical persistent pollutants are smaller (<=700 Da)"
    
    # Define recursive SMARTS for the characteristic scaffolds.
    # The patterns are designed to capture the main aromatic core even when substituted.
    biphenyl_pattern = Chem.MolFromSmarts("[$(c1ccccc1)-$(c2ccccc2)]")
    dibenzodioxin_pattern = Chem.MolFromSmarts("[$(c1ccc2Oc3ccccc3O2c1)]")
    dibenzofuran_pattern = Chem.MolFromSmarts("[$(c1ccc2Oc3ccccc3c2c1)]")
    
    scaffold_matches = []
    if mol.HasSubstructMatch(biphenyl_pattern):
        scaffold_matches.append("biphenyl")
    if mol.HasSubstructMatch(dibenzodioxin_pattern):
        scaffold_matches.append("dibenzodioxin")
    if mol.HasSubstructMatch(dibenzofuran_pattern):
        scaffold_matches.append("dibenzofuran")
        
    if not scaffold_matches:
        return False, "No recognized aromatic scaffold (biphenyl, dibenzodioxin, or dibenzofuran) found"
    
    # Build a reason string summarizing the molecule's features.
    reason = (
        f"Molecule has a {' and '.join(scaffold_matches)} scaffold with {len(halogen_atoms)} halogen substituents, "
        f"a molecular weight of {mol_wt:.1f} Da"
    )
    
    return True, reason

# Example usage (for testing):
if __name__ == '__main__':
    # Test with an example: 4,4'-dichlorobiphenyl
    test_smiles = "Clc1ccc(cc1)-c1ccc(Cl)cc1"
    result, msg = is_polychlorinated_dibenzodioxines_and_related_compounds(test_smiles)
    print(result, msg)