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

    # Check for specific tertiary amine patterns
    tertiary_amine_smarts = ["[N;H0;X3]([C])[C]", "[N;H0;X3]([C])([C])[C]", "[N;H0;X3]([C])([c])[C]", "[N;H0;X3]([C])([c])([c])"]
    for smarts in tertiary_amine_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a tertiary nitrogen atom with hydrocarbyl groups attached"

    return True, "Matches tertiary amine criteria based on structural analysis"