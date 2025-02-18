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
    A dicarboxylic acid has exactly two carboxy (-COOH) groups that are not part of an anhydride.

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

    # Define carboxy group pattern (accounts for -COOH or -COO-)
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

    # Find all matches
    matches = mol.GetSubstructMatches(carboxy_pattern)
    
    # Check for exactly two carboxy groups
    if len(matches) != 2:
        return False, f"Found {len(matches)} carboxy groups, need exactly 2"

    # Collect all atoms in the carboxy groups
    atoms_in_groups = set()
    for match in matches:
        atoms_in_groups.update(match)  # match is a tuple of atom indices

    # Get all oxygen atoms in the molecule
    o_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    
    # Count oxygen atoms that are part of the carboxy groups
    shared_os = sum(1 for o in o_atoms if o in atoms_in_groups)
    
    # Valid dicarboxylic acid must have exactly 4 oxygen atoms in the carboxy groups (no sharing)
    if shared_os == 4:
        return True, "Contains two distinct carboxy groups"
    else:
        return False, "Carboxy groups share oxygen atoms (likely anhydride or conjugated)"