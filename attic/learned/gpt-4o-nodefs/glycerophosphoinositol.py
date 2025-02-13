"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for common inositol ring variations
    inositol_patterns = [
        Chem.MolFromSmarts("C1([C@H]([C@@H]([C@@H]([C@H]([C@H]1O)O)O)O)O)O"),  # inositol with varied stereochemistry
        Chem.MolFromSmarts("C1(C(C(C(C(C1O)O)O)O)O)O")  # common inositol pattern
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in inositol_patterns):
        return False, "No matching inositol ring structure found"

    # Check for phosphate group attached in a phosphatidyl form
    phosphate_patterns = [
        Chem.MolFromSmarts("O[P@](=O)(O)O"),  # specific configuration phosphate
        Chem.MolFromSmarts("OP(=O)(O)O")  # general phosphate group
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns):
        return False, "No phosphate group found"

    # Check for glycerol linkage indicative of lipid backbones
    glycerol_patterns = [
        Chem.MolFromSmarts("[C@H](CO)CO"),  # basic glycerol structure
        Chem.MolFromSmarts("C(O)C(O)CO")  # varied configuration
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in glycerol_patterns):
        return False, "No glycerol backbone found"

    # Identify presence of hydrophobic fatty chains
    chain_patterns = [
        Chem.MolFromSmarts("C[C@H](C)CC"), # looking for generic aliphatic chains
        Chem.MolFromSmarts("CCCCCCCC")  # presence of a carbon chain
    ]
    chain_matches = sum(mol.HasSubstructMatch(pattern) for pattern in chain_patterns)
    if chain_matches < 1:
        return False, "Inadequate hydrophobic fatty acid chains detected"

    return True, "Contains essential glycerophosphoinositol structures with inositol ring, phosphate group, and probable fatty acid chains"

# Example usage:
# result, message = is_glycerophosphoinositol("P(O[C@H]1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCCC/C=C\C/C=C\C/C=C\CCCCC)COC(=O)CCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)(O)=O")
# print(result, message)