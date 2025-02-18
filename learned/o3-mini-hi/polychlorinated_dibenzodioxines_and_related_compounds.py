"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Organochlorine (and/or organobromine) compounds that are polychlorinated dibenzodioxines and related compounds.
Definition:
  - They are persistent organic pollutants containing Cl and/or Br.
  - They must contain a minimum number of aromatic halogen substituents (halogen directly attached to an aromatic carbon):
      • biphenyl or diaryl ether scaffold: at least 2 aromatic halogens.
      • dibenzofuran scaffold: at least 3 aromatic halogens.
      • dibenzodioxin scaffold: at least 4 aromatic halogens.
  - Their molecular weight must be no more than 700 Da.
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
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxines and related compounds.

    Criteria:
      1. The molecule must contain a minimum number of aromatic halogen substituents (Cl (Z=17) or Br (Z=35) attached to an aromatic carbon).
         The minimum required depends on the aromatic scaffold that is present:
           - biphenyl or diaryl ether: require >= 2 aromatic halogens.
           - dibenzofuran: require >= 3 aromatic halogens.
           - dibenzodioxin: require >= 4 aromatic halogens.
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
    # Convert SMILES to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First check: molecular weight must be <=700 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 700:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da); should be <= 700 Da"
    
    # Count only halogen atoms (Cl or Br) that are attached to an aromatic carbon.
    aromatic_halogen_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in (17, 35):  # Cl or Br
            # Check if at least one neighbor is aromatic.
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic():
                    aromatic_halogen_count += 1
                    break

    if aromatic_halogen_count < 2:
        return False, "Insufficient aromaticly bound halogen substituents (need at least 2)"

    # Define scaffold SMARTS patterns.
    scaffold_smarts = {
        "biphenyl": "c1ccccc1-c2ccccc2",
        "diaryl ether": "c1ccccc1Oc2ccccc2",
        "dibenzodioxin": "c1ccc2Oc3ccccc3O2c1",
        "dibenzofuran": "c1ccc2oc3ccccc3c2c1"
    }
    # Minimum aromatic halogen counts required by the scaffold type.
    scaffold_halogen_requirements = {
        "biphenyl": 2,
        "diaryl ether": 2,
        "dibenzodioxin": 4,
        "dibenzofuran": 3
    }
    
    scaffolds_found = []
    # Loop over scaffolds and check for presence.
    for scaf_name, smarts in scaffold_smarts.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # Should not occur.
        if mol.HasSubstructMatch(pattern):
            scaffolds_found.append(scaf_name)
    
    if not scaffolds_found:
        return False, "No recognized aromatic scaffold (biphenyl, diaryl ether, dibenzodioxin, or dibenzofuran) found"
    
    # Determine the most “demanding” (highest) halogen requirement among the scaffolds detected.
    required_halogen = 0
    for scaf in scaffolds_found:
        req = scaffold_halogen_requirements.get(scaf, 2)
        if req > required_halogen:
            required_halogen = req

    if aromatic_halogen_count < required_halogen:
        return False, (
            f"Scaffold(s) {', '.join(scaffolds_found)} detected but only {aromatic_halogen_count} aromatic halogen(s) found; "
            f"need at least {required_halogen} for this scaffold type"
        )
    
    reason = (
        f"Molecule has a {', '.join(scaffolds_found)} scaffold, "
        f"{aromatic_halogen_count} aromatic halogen(s), and a molecular weight of {mol_wt:.1f} Da"
    )
    return True, reason

# Example usage for testing.
if __name__ == '__main__':
    # Test with one known true-positive: 4,4'-dichlorobiphenyl
    test_smiles = "Clc1ccc(cc1)-c1ccc(Cl)cc1"
    result, msg = is_polychlorinated_dibenzodioxines_and_related_compounds(test_smiles)
    print(result, msg)