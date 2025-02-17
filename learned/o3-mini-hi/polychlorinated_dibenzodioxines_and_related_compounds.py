"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Organochlorine compounds that are polychlorinated dibenzodioxines and related compounds.
Definition:
  - They are organochlorine (and/or organobromine) compounds that are persistent organic pollutants.
  - They must contain a minimum number of halogen atoms (Cl or Br) – but the minimum depends
    on their aromatic scaffold:
      • For a biphenyl or diaryl ether scaffold: at least 2 halogens.
      • For a dibenzofuran scaffold: at least 3 halogens.
      • For a dibenzodioxin scaffold: at least 4 halogens.
  - Their molecular weight should be no more than 700 Da.
  - They must contain (at least one) of the following aromatic scaffolds:
         • Biphenyl (two benzene rings connected by a single bond)
         • Diaryl ether (two benzene rings connected by an oxygen atom)
         • Dibenzodioxin (two benzene rings connected by two oxygens in a dioxin motif)
         • Dibenzofuran (two benzene rings connected by an oxygen in a furan motif)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated
    dibenzodioxines and related compounds.

    Criteria:
      1. The molecule must contain a minimum number of halogen atoms (Cl (Z=17) or Br (Z=35)).
         However, the minimum required depends on the scaffold that is present:
           - biphenyl or diaryl ether: require at least 2 halogens.
           - dibenzofuran: require at least 3 halogens.
           - dibenzodioxin: require at least 4 halogens.
      2. The molecular weight must be <= 700 Da.
      3. The molecule must contain at least one of the characteristic aromatic scaffolds:
             • Biphenyl: SMARTS "c1ccccc1-c2ccccc2"
             • Diaryl ether: SMARTS "c1ccccc1Oc2ccccc2"
             • Dibenzodioxin: SMARTS "c1ccc2Oc3ccccc3O2c1"
             • Dibenzofuran: SMARTS "c1ccc2oc3ccccc3c2c1"

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets the combined criteria, False otherwise.
        str: Explanation for the classification decision.
    """

    # Convert SMILES to RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count halogen substituents (only Cl (17) or Br (35))
    halogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in (17, 35)]
    total_halogens = len(halogen_atoms)
    if total_halogens < 2:
        return False, "Insufficient halogen atoms (need at least 2 Cl or Br substituents)"
    
    # Verify molecular weight <= 700 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 700:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da); should be <= 700 Da"
    
    # Define scaffold SMARTS patterns.
    scaffold_smarts = {
        "biphenyl": "c1ccccc1-c2ccccc2",
        "diaryl ether": "c1ccccc1Oc2ccccc2",
        "dibenzodioxin": "c1ccc2Oc3ccccc3O2c1",
        "dibenzofuran": "c1ccc2oc3ccccc3c2c1"
    }
    # Minimum halogen counts required by the scaffold type:
    scaffold_halogen_requirements = {
        "biphenyl": 2,
        "diaryl ether": 2,
        "dibenzodioxin": 4,
        "dibenzofuran": 3
    }
    
    scaffolds_found = []
    valid_scaffold = False
    # Loop over scaffolds and check for presence.
    for scaf, smarts in scaffold_smarts.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # Should not happen.
        if mol.HasSubstructMatch(pattern):
            scaffolds_found.append(scaf)
    
    if not scaffolds_found:
        return False, "No recognized aromatic scaffold (biphenyl, diaryl ether, dibenzodioxin, or dibenzofuran) found"
    
    # Now, determine the (highest) halogen requirement implied by the matched scaffolds.
    # (If a molecule matches more than one scaffold, we use the most “demanding” requirement.)
    required_halogens = 0
    for scaf in scaffolds_found:
        req = scaffold_halogen_requirements.get(scaf, 2)
        if req > required_halogens:
            required_halogens = req

    if total_halogens < required_halogens:
        return False, (
            f"Scaffold(s) {', '.join(scaffolds_found)} detected but only {total_halogens} halogen(s) present; "
            f"need at least {required_halogens} for this scaffold type"
        )
    
    # If all criteria are met, build a summary reason.
    reason = (
        f"Molecule has a {', '.join(scaffolds_found)} scaffold, "
        f"{total_halogens} halogen(s), and a molecular weight of {mol_wt:.1f} Da"
    )
    return True, reason

# Example usage for testing.
if __name__ == '__main__':
    # Test with one known true-positive: 4,4'-dichlorobiphenyl
    test_smiles = "Clc1ccc(cc1)-c1ccc(Cl)cc1"
    result, msg = is_polychlorinated_dibenzodioxines_and_related_compounds(test_smiles)
    print(result, msg)