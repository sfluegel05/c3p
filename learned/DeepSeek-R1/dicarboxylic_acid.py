"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:17895 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid has exactly two carboxy (-COOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxy group pattern (accounts for different tautomers/protonation states)
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")  # Matches -COOH or -COO-

    # Find all matches
    matches = mol.GetSubstructMatches(carboxy_pattern)
    
    # Check for exactly two carboxy groups
    if len(matches) != 2:
        return False, f"Found {len(matches)} carboxy groups, need exactly 2"

    # Ensure the two groups are not conjugated (e.g., not part of an anhydride)
    # Check if the two carboxylic acids share oxygen atoms (common in anhydrides)
    atoms_in_groups = set()
    for match in matches:
        atoms_in_groups.update(match.GetIndices())
    
    # If any oxygen is shared between the two groups, it's likely an anhydride
    o_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    shared_os = len([o for o in o_atoms if o in atoms_in_groups])
    if shared_os < 4:  # Each COOH has 2 O atoms, total 4 if not shared
        return True, "Contains two distinct carboxy groups"
    else:
        return False, "Carboxy groups share oxygen atoms (likely anhydride)"