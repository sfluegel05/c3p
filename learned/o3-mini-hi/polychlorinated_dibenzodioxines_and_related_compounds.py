"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Organochlorine compounds that are polychlorinated dibenzodioxines and structurally related entities.
These include polychlorinated dibenzofurans as well as polychlorinated and polybrominated biphenyls.
A molecule must have at least two chlorine or bromine atoms and contain one of the characteristic aromatic scaffolds.
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxines and related compounds.
    
    This classification is based on two conditions:
      1. The molecule must contain aromatic halogen substituents (chlorine or bromine) in sufficient number.
      2. The molecule must contain at least one of the following structural scaffolds:
           a. Biphenyl: two benzene rings connected by a single bond.
           b. Dibenzodioxin: two benzene rings linked via two oxygen bridges.
           c. Dibenzofuran: two benzene rings sharing an oxygen in a fused ring system.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule matches criteria, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Count halogens (chlorine atomic num 17 and bromine atomic num 35)
    halogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in (17, 35)]
    if len(halogen_atoms) < 2:
        return False, "Insufficient halogen atoms (need at least 2 Cl or Br substituents)"
        
    # Define SMARTS patterns for different scaffolds
    # Biphenyl scaffold: two benzene rings connected by a single bond.
    biphenyl_smarts = "c1ccccc1-c2ccccc2"
    biphenyl_pattern = Chem.MolFromSmarts(biphenyl_smarts)
    
    # Dibenzodioxin scaffold: two benzene rings connected via two oxygen bridges.
    dibenzodioxin_smarts = "c1cc2Oc3ccccc3O2c1"
    dibenzodioxin_pattern = Chem.MolFromSmarts(dibenzodioxin_smarts)
    
    # Dibenzofuran scaffold: two benzene rings fused through an oxygen.
    dibenzofuran_smarts = "c1ccc2Oc3ccccc3c2c1"
    dibenzofuran_pattern = Chem.MolFromSmarts(dibenzofuran_smarts)
    
    # Check for at least one scaffold match
    scaffold_found = False
    scaffold_matches = []
    if mol.HasSubstructMatch(biphenyl_pattern):
        scaffold_found = True
        scaffold_matches.append("biphenyl")
    if mol.HasSubstructMatch(dibenzodioxin_pattern):
        scaffold_found = True
        scaffold_matches.append("dibenzodioxin")
    if mol.HasSubstructMatch(dibenzofuran_pattern):
        scaffold_found = True
        scaffold_matches.append("dibenzofuran")
        
    if not scaffold_found:
        return False, "No recognized aromatic scaffold (biphenyl, dibenzodioxin or dibenzofuran) found"
    
    # If we reach here, the compound has enough halogen substituents and a matching scaffold.
    reason = "Molecule has a %s scaffold with %d halogen substituents" % (
        " and ".join(scaffold_matches), len(halogen_atoms))
    return True, reason

# Example usage (for testing):
if __name__ == '__main__':
    # Test with 4,4'-dichlorobiphenyl example
    test_smiles = "Clc1ccc(cc1)-c1ccc(Cl)cc1"
    result, msg = is_polychlorinated_dibenzodioxines_and_related_compounds(test_smiles)
    print(result, msg)