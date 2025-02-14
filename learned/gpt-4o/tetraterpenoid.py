"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (35 <= c_count <= 45):
        return False, f"Expected around 40 carbon atoms, found {c_count}"

    # Check for conjugated double bonds pattern
    # SMARTS for conjugated polyene segment (repeated C=C-C-C); flexible number of repetitions
    polyene_pattern = Chem.MolFromSmarts("C=C(-C)=C")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No conjugated polyene backbone found"

    # Check for presence of typical functional groups (e.g., -OH, =O)
    functional_patterns = [
        Chem.MolFromSmarts("O"),  # Hydroxyl
        Chem.MolFromSmarts("C=O"),  # Carbonyl
        Chem.MolFromSmarts("O[C@H]")  # Epoxide (stereochemistry typical for carotenoids)
    ]
    found_functional_groups = any(mol.HasSubstructMatch(pat) for pat in functional_patterns)
    if not found_functional_groups:
        return False, "No common functional groups (e.g., hydroxyl, carbonyl) found"

    return True, "Structure consistent with a tetraterpenoid (typical C40 backbone with isoprene units and carotenoid-like functional groups)"