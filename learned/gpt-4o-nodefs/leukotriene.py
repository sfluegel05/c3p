"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    Leukotrienes have polyene structures, usually derived from arachidonic acid,
    containing conjugated double bonds, carboxyl groups, and may include hydroxyl
    and sulfur-containing groups.

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

    # Revise carbon count to encompass broader variations in chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 58):
        return False, f"Unusual carbon count: found {c_count} carbons"

    # Conjugated double bond patterns - expanded possibilities
    conjugated_patterns = [
        Chem.MolFromSmarts("[#6]=[#6]-[#6]=[#6]-[#6]=[#6]"),
        Chem.MolFromSmarts("[#6]-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]"),
        Chem.MolFromSmarts("[#6]=[#6]-[#6]=[#6]"),
        Chem.MolFromSmarts("[#6]-[#6]=[#6]-[#6]=[#6]"),
        Chem.MolFromSmarts("[#6]-[#6]=[#6]-[#6]-[#6]=[#6]"),
        Chem.MolFromSmarts("[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]"),
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in conjugated_patterns):
        return False, "No suitable conjugated double bond pattern found"

    # Carboxyl groups presence
    carboxyl_patterns = [Chem.MolFromSmarts("C(=O)[OH]"), Chem.MolFromSmarts("C(=O)[O-]")]
    if not any(mol.HasSubstructMatch(pattern) for pattern in carboxyl_patterns):
        return False, "No carboxyl group (COOH or COO-) detected"

    # Check for hydroxyl groups (flexible count)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if len(mol.GetSubstructMatches(hydroxyl_pattern)) == 0:
        return False, "No hydroxyl group (OH) found"

    # Consider sulfur - optional due to variety in leukotriene types
    sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    # Optional validation if there are distinguishing sulfur patterns

    return True, "Structure matches typical features of leukotrienes"