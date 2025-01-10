"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    Leukotrienes typically have polyene structures derived from arachidonic acid 
    and may include hydroxyl, carboxyl, and sulfur-containing groups.

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

    # General carbon count flexibility - updated range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 32):
        return False, f"Unusual carbon count: found {c_count} carbons"

    # Conjugated double bond patterns - extending to broader possibilities
    conjugated_patterns = [
        Chem.MolFromSmarts("[#6]=[#6]-[#6]=[#6]-[#6]=[#6]"),
        Chem.MolFromSmarts("[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]"),
        Chem.MolFromSmarts("[#6]=[#6]-[#6]-[#6]=[#6]-[#6]=[#6]"),
        Chem.MolFromSmarts("[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]")  # Account for terminal alkene groups
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in conjugated_patterns):
        return False, "No suitable conjugated double bond pattern found"

    # At least one carboxyl group
    carboxyl_patterns = [Chem.MolFromSmarts("C(=O)[OH]"), Chem.MolFromSmarts("C(=O)[O-]")]
    if not any(mol.HasSubstructMatch(pattern) for pattern in carboxyl_patterns):
        return False, "No carboxyl group (COOH or COO-) detected"

    # Hydroxyl group qualification
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    alcohol_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if alcohol_count == 0:
        return False, "No hydroxyl group (OH) found"

    # Qualify based on sulfur content, but allow for absence as well
    sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if sulfur_count > 0:
        return True, "Classes as cysteinyl leukotriene with sulfur present"
    
    return True, "Structure matches typical features of leukotrienes"