"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:35617 tertiary amine
A compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define tertiary amine pattern
    tertiary_amine_pattern = Chem.MolFromSmarts("[N;H0;X3]")

    # Check for tertiary nitrogen atoms
    tertiary_n_atoms = mol.GetSubstructMatches(tertiary_amine_pattern)
    if not tertiary_n_atoms:
        return False, "No tertiary nitrogen atoms found"

    # Check for hydrocarbyl groups (alkyl or aryl) attached to tertiary nitrogen atoms
    hydrocarbyl_pattern = Chem.MolFromSmarts("[N;H0;X3][C]")
    has_hydrocarbyl = any(mol.HasSubstructMatch(hydrocarbyl_pattern))
    if not has_hydrocarbyl:
        return False, "No hydrocarbyl groups attached to tertiary nitrogen atoms"

    # Count number of hydrocarbyl groups attached to each tertiary nitrogen
    num_hydrocarbyl_groups = [len(mol.GetAtomWithIdx(atom_idx).GetNeighbors()) for atom_idx in tertiary_n_atoms]
    if any(n != 3 for n in num_hydrocarbyl_groups):
        return False, "Not all tertiary nitrogen atoms have exactly three hydrocarbyl groups attached"

    # Check if the hydrocarbyl groups are valid (alkyl or aryl)
    valid_hydrocarbyl_pattern = Chem.MolFromSmarts("[C;H3,H2,H1,H0]~[C;H3,H2,H1,H0]")
    for atom_idx in tertiary_n_atoms:
        tertiary_n = mol.GetAtomWithIdx(atom_idx)
        neighbors = tertiary_n.GetNeighbors()
        for neighbor in neighbors:
            if not mol.HasSubstructMatch(valid_hydrocarbyl_pattern, atomIdxList=[neighbor.GetIdx()]):
                return False, "At least one hydrocarbyl group attached to a tertiary nitrogen is not alkyl or aryl"

    return True, "Contains a tertiary nitrogen atom with three hydrocarbyl groups attached"