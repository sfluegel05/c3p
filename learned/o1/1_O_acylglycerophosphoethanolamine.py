"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    This class of molecules has a glycerol backbone with:
    - An acyl chain attached via an ester bond at position 1
    - A free hydroxyl group at position 2
    - A phosphoethanolamine group attached at position 3
    - No other acyl chains or substitutions on the glycerol backbone

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    from rdkit import Chem

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphoethanolamine group
    phosphoethanolamine_smarts = 'O[P](=O)(O)OCC[NH2,NH3+]'
    phosphoethanolamine_pattern = Chem.MolFromSmarts(phosphoethanolamine_smarts)
    if phosphoethanolamine_pattern is None:
        return False, "Error in phosphoethanolamine SMARTS pattern"

    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found"

    # Check for glycerol backbone with acyl chain at position 1 and free hydroxyl at position 2
    glycerol_smarts = '[C@@H]([O][CH2][P](=O)([O])[O][CH2][CH2][NH2,NH3+])[CH](O)[CH2][O][C](=O)[C]'
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    if glycerol_pattern is None:
        return False, "Error in glycerol backbone SMARTS pattern"

    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with correct substitutions not found"

    # Ensure there is only one ester group (acyl chain)
    ester_smarts = 'C(=O)O[C][CH]'
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Error in ester group SMARTS pattern"

    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, expected exactly 1"

    # Check for free hydroxyl group at position 2
    free_oh_smarts = '[C@H](O)[CH2][O][P](=O)(O)[O][CH2][CH2][NH2,NH3+]'
    free_oh_pattern = Chem.MolFromSmarts(free_oh_smarts)
    if free_oh_pattern is None:
        return False, "Error in free hydroxyl SMARTS pattern"

    if not mol.HasSubstructMatch(free_oh_pattern):
        return False, "No free hydroxyl group at position 2 found"

    # Ensure there are no additional acyl chains or substitutions
    # Count all ester bonds in the molecule
    ester_bond_pattern = Chem.MolFromSmarts('C(=O)O[C]')
    ester_bond_matches = mol.GetSubstructMatches(ester_bond_pattern)
    if len(ester_bond_matches) != 1:
        return False, "Additional ester bonds found, molecule may have extra acyl chains"

    return True, "Molecule is a 1-O-acylglycerophosphoethanolamine"