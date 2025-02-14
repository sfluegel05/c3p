"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: CHEBI:26189 polymer
"""

from rdkit import Chem

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    This is done by analyzing atom environments to detect repeating units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a polymer, False otherwise
        str: Reason for classification
    """
    from collections import defaultdict

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Dictionary to count substructures
    substruct_counts = defaultdict(int)
    radius = 2  # radius for neighborhoods

    for atom_idx in range(mol.GetNumAtoms()):
        # Get the atom environment for the given radius
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)
        if not env:
            continue  # Skip if no environment is found (e.g., for isolated atoms)
        amap = {}
        submol = Chem.PathToSubmol(mol, env, atomMap=amap)
        # Generate the canonical SMILES of the substructure
        smiles_submol = Chem.MolToSmiles(submol, canonical=True)
        substruct_counts[smiles_submol] += 1

    # Find the maximum count of any substructure
    if not substruct_counts:
        return False, "No substructures found for analysis"

    max_count = max(substruct_counts.values())

    # Threshold for repeating units
    repeating_unit_threshold = 5  # This value may need adjustment

    if max_count >= repeating_unit_threshold:
        return True, f"Molecule has repeating substructures occurring {max_count} times, likely a polymer"
    else:
        return False, "No significant repeating units found"