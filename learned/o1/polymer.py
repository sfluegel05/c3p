"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: CHEBI:26189 polymer
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is characterized by high molecular weight and large size.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a polymer, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate exact molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Get the number of atoms
    num_atoms = mol.GetNumAtoms()

    # Define thresholds for polymers (these values are somewhat arbitrary)
    molecular_weight_threshold = 1000  # Da
    atom_count_threshold = 100

    # Check if molecular weight is above the threshold
    if mol_wt < molecular_weight_threshold:
        return False, f"Molecular weight {mol_wt:.2f} Da is too low for a polymer"

    # Check if number of atoms is above the threshold
    if num_atoms < atom_count_threshold:
        return False, f"Molecule has only {num_atoms} atoms, too small for a polymer"

    # Since detecting repeating units is complex, and cannot be accurately done,
    # we will assume that molecules with very high molecular weight and large number of atoms
    # are likely polymers

    return True, f"Molecule has high molecular weight ({mol_wt:.2f} Da) and large size ({num_atoms} atoms), likely a polymer"