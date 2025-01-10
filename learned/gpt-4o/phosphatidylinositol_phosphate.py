"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Identifies phosphatidylinositol phosphates, which are a class of phospholipids containing a glycerol backbone,
    an inositol ring, and one or more phosphate groups attached to the inositol ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: Identification result, True if molecule is a phosphatidylinositol phosphate
        str: Explanation for the classification
    """
    
    # Create rdkit Molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone pattern: matching glycerol with ester bonds
    glycerol_pattern = Chem.MolFromSmarts("C(COC(=O))COC(=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with ester bonds found"

    # Inositol ring pattern: a cyclohexane ring with hydroxyl groups
    inositol_pattern = Chem.MolFromSmarts("C1(CO)C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Phosphate group pattern: looking for one or more phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("O[P](=O)(O)O")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_count < 1:
        return False, "Less than one phosphate group found"

    return True, "Contains glycerol backbone with inositol ring and one or more phosphate groups"

# Example test
result, reason = is_phosphatidylinositol_phosphate("CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC\C=C/CCCCCCCC")
print(result, reason)  # Expected output: True, with detailed reason