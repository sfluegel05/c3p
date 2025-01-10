"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    Leukotrienes feature conjugated polyene structures typically derived from
    arachidonic acid and are often associated with hydroxyl or carboxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexible carbon count (around 20 is typical but not rigid)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 26):
        return False, f"Unusual carbon count: found {c_count} carbons"

    # Detect various configurations of conjugated double bonds
    conjugated_patterns = [
        Chem.MolFromSmarts("C=C-C=C-C=C"),
        Chem.MolFromSmarts("C=C-C=C-C=C-C=C"),
        Chem.MolFromSmarts("C=C-C-C=C-C=C")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in conjugated_patterns):
        return False, "No suitable conjugated double bond pattern found"

    # Check for at least one carboxyl group (COOH) or carboxylate ion (COO-)
    carboxyl_patterns = [Chem.MolFromSmarts("C(=O)O"), Chem.MolFromSmarts("C(=O)[O-]")]
    if not any(mol.HasSubstructMatch(pattern) for pattern in carboxyl_patterns):
        return False, "No carboxyl group (COOH or COO-) detected"

    # Allow for multiple hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    alcohol_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if alcohol_count == 0:
        return False, "No hydroxyl group (OH) found"

    # Presence of sulfur is optional for cysteinyl leukotrienes like LTC4
    sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if sulfur_count > 0:
        return True, "Classes as cysteinyl leukotriene with sulfur present"

    return True, "Structure matches typical features of leukotrienes"