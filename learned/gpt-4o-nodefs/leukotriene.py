"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    Leukotrienes are characterized by a conjugated polyene structure, typically
    derived from arachidonic acid, with some hydroxyl or carboxyl groups.

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

    # Count carbon atoms (Expecting 20 for leukotrienes)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Carbon count mismatch: found {c_count} carbons, expected 20"

    # Detect conjugated double bonds (at least 3 consecutive double bonds)
    # Simplified assumption here: we are looking for a common pattern of conjugated dienes
    conjugated_diene_pattern = Chem.MolFromSmarts("C=C-C=C")
    conjugated_diene_matches = mol.GetSubstructMatches(conjugated_diene_pattern)
    if len(conjugated_diene_matches) < 3:
        return False, "Lacking sufficient conjugated double bond structures"

    # Check for a carboxyl group (COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group (COOH) detected"

    # (Optional) Check for alcohol groups: -OH
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No hydroxyl group (OH) detected"

    # Check for sulfur atom in cases of cysteinyl leukotrienes like LTC4
    sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if sulfur_count > 0:
        return True, "Classes as cysteinyl leukotriene with sulfur present"

    return True, "Structure matches typical features of leukotrienes"