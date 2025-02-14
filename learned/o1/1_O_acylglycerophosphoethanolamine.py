"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    This class of molecules has a glycerol backbone with:
    - An acyl chain attached via an ester bond at position 1
    - A free hydroxyl group at position 2
    - A phosphoethanolamine group at position 3

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

    # Define SMARTS patterns
    # Phosphoethanolamine group
    phosphoethanolamine_pattern = Chem.MolFromSmarts('COP(O)(=O)OCCN')
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found"

    # Glycerol backbone with acyl group at position 1 and free OH at position 2
    # Acyl group at position 1
    acyl_glycerol_pattern = Chem.MolFromSmarts('[$([C@@H](CO[P](=O)(O)OCCN)OC(=O)[#6])][$([C@H](O)[O])]')  
    if not mol.HasSubstructMatch(acyl_glycerol_pattern):
        return False, "No glycerol backbone with acyl group at position 1 and free OH at position 2 found"

    # Check for only one acyl chain attached via ester bond
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C;H]')  # Ester linkage to primary carbon
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester linkage, found {len(ester_matches)}"

    # Check for free hydroxyl group at position 2 of glycerol
    free_oh_pattern = Chem.MolFromSmarts('[C@H](O)[CH2]O[P](=O)(O)OCCN')
    if not mol.HasSubstructMatch(free_oh_pattern):
        return False, "No free hydroxyl group at position 2 of glycerol found"

    return True, "Molecule is a 1-O-acylglycerophosphoethanolamine"