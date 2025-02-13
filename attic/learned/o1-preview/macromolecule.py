"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is defined as a molecule of high relative molecular mass,
    comprising multiple repetitions of units derived from small molecules.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight is {mol_wt:.2f} Da, which is below the macromolecule threshold"

    # Count potential repeating units
    # Common monomer units: amino acids, monosaccharides, nucleotides
    # Define SMARTS patterns for peptide bonds, glycosidic bonds, and phosphodiester bonds

    # Peptide bond pattern (amide linkage between amino acids)
    peptide_bond = Chem.MolFromSmarts("C(=O)N")
    num_peptide_bonds = len(mol.GetSubstructMatches(peptide_bond))

    # Glycosidic bond pattern (ether linkage between sugars)
    glycosidic_bond = Chem.MolFromSmarts("[C;H1,H2,H3][O][C;H1,H2,H3]")
    num_glycosidic_bonds = len(mol.GetSubstructMatches(glycosidic_bond))

    # Phosphodiester bond pattern (linkage in nucleic acids)
    phosphodiester_bond = Chem.MolFromSmarts("P(=O)(O)O[C;H1,H2,H3]")
    num_phosphodiester_bonds = len(mol.GetSubstructMatches(phosphodiester_bond))

    # Check if the molecule has multiple repeating units
    total_repeating_bonds = num_peptide_bonds + num_glycosidic_bonds + num_phosphodiester_bonds

    if total_repeating_bonds >= 5:
        return True, f"Molecule has high molecular weight and contains {total_repeating_bonds} repeating linkage units characteristic of macromolecules"

    # Alternatively, check for large number of similar substructures
    # Use fingerprints to detect repeating patterns
    fingerprints = {}
    for atom in mol.GetAtoms():
        # Get local environment (e.g., atom and neighbors)
        env = Chem.PathToSubmol(mol, [atom.GetIdx()])
        smi = Chem.MolToSmiles(env, canonical=True)
        fingerprints[smi] = fingerprints.get(smi, 0) + 1

    # If any substructure appears multiple times, it may indicate repeating units
    max_repeat = max(fingerprints.values())
    if max_repeat >= 5:
        return True, f"Molecule has high molecular weight and contains repeating units detected by substructure frequency"

    return False, "Molecule does not meet the criteria for a macromolecule"

# __metadata__ section can remain the same as before