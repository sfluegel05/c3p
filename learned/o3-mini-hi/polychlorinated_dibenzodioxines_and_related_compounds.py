"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Organochlorine compounds that are polychlorinated dibenzodioxines and related compounds.
These include polychlorinated dibenzofurans and polychlorinated (or polybrominated) biphenyls.
The criteria are:
  - The molecule must contain at least 2 chlorine or bromine atoms.
  - It must contain one of the characteristic aromatic scaffolds –
       (a) biphenyl,
       (b) dibenzodioxin,
       (c) dibenzofuran.
  - In addition, it should be of similar size to typical persistent organic pollutants.
    Here we require the molecular weight to be below 600 Da and that the molecule has exactly 2 aromatic rings.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxines and related compounds.
    The classification combines three sets of criteria:
      1. The molecule must have at least 2 halogen atoms (Cl or Br).
      2. The underlying aromatic scaffold must be one of:
             biphenyl, dibenzodioxin, or dibenzofuran.
         We use recursive SMARTS patterns to allow extra substituents.
      3. The molecule should be similar in size to persistent organic pollutants:
             (a) Molecular weight below 600 Da.
             (b) Exactly 2 aromatic rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
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
    if mol_wt >= 600:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da); typical pollutants are smaller (<600 Da)"
    
    # Count the number of aromatic rings.
    aromatic_ring_count = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_ring_count != 2:
        return False, f"Unexpected number of aromatic rings ({aromatic_ring_count}); expected exactly 2"
    
    # Define recursive SMARTS for the characteristic scaffolds.
    # Using the $() syntax allows for extra substituents on the aromatic rings.
    biphenyl_smarts = "[$(c1ccccc1)-$(c2ccccc2)]"
    dibenzodioxin_smarts = "[$(c1ccc2Oc3ccccc3O2c1)]"
    dibenzofuran_smarts = "[$(c1ccc2Oc3ccccc3c2c1)]"
    
    biphenyl_pattern = Chem.MolFromSmarts(biphenyl_smarts)
    dibenzodioxin_pattern = Chem.MolFromSmarts(dibenzodioxin_smarts)
    dibenzofuran_pattern = Chem.MolFromSmarts(dibenzofuran_smarts)
    
    scaffold_matches = []
    if mol.HasSubstructMatch(biphenyl_pattern):
        scaffold_matches.append("biphenyl")
    if mol.HasSubstructMatch(dibenzodioxin_pattern):
        scaffold_matches.append("dibenzodioxin")
    if mol.HasSubstructMatch(dibenzofuran_pattern):
        scaffold_matches.append("dibenzofuran")
        
    if not scaffold_matches:
        return False, "No recognized aromatic scaffold (biphenyl, dibenzodioxin or dibenzofuran) found"
    
    # Build a reason string summarizing the molecule's features.
    reason = ("Molecule has a %s scaffold with %d halogen substituents, "
              "molecular weight of %.1f Da and %d aromatic rings") % (
                  " and ".join(scaffold_matches), len(halogen_atoms), mol_wt, aromatic_ring_count)
    
    return True, reason

# Example usage (for testing):
if __name__ == '__main__':
    # Test with one example: 4,4'-dichlorobiphenyl
    test_smiles = "Clc1ccc(cc1)-c1ccc(Cl)cc1"
    result, msg = is_polychlorinated_dibenzodioxines_and_related_compounds(test_smiles)
    print(result, msg)