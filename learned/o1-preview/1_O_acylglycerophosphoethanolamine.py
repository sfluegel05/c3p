"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.

    A 1-O-acylglycerophosphoethanolamine is a glycerophosphoethanolamine having an O-acyl substituent at the 1-position of the glycerol fragment.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the 1-O-acylglycerophosphoethanolamine backbone
    pattern_smarts = """
    [C@@H]([O][C](=O)[C])[C@H](O)COP(=O)(O)OCCN
    """
    pattern = Chem.MolFromSmarts(pattern_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    if not mol.HasSubstructMatch(pattern):
        return False, "Does not match 1-O-acylglycerophosphoethanolamine substructure"

    # Ensure there is exactly one acyl chain attached via ester linkage
    ester_pattern = Chem.MolFromSmarts("[C@@H](O)COC(=O)[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected one acyl chain attached via ester linkage, found {len(ester_matches)}"

    # Check for free hydroxyl group at position 2 of glycerol backbone
    hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)COP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No free hydroxyl group at position 2 of glycerol backbone"

    # Ensure there are no additional acyl chains or ether linkages
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    acyl_chains = len(mol.GetSubstructMatches(acyl_chain_pattern))
    if acyl_chains != 1:
        return False, f"Expected one acyl chain, found {acyl_chains}"

    # Check for phosphoethanolamine group specifically
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "Phosphoethanolamine group not found at position 3"

    return True, "Matches 1-O-acylglycerophosphoethanolamine structure"