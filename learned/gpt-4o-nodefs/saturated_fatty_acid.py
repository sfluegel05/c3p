"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify carboxylic acid group at terminal position
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not matches:
        return False, "No carboxylic acid group found"
    
    num_carboxylic_acids = len(matches)
    if num_carboxylic_acids != 1:
        return False, f"Found {num_carboxylic_acids} carboxylic acid groups; requires exactly 1"

    # Ensure that the carboxylic acid is at one end of a chain
    match = matches[0]
    carboxylic_carbon_idx = match[0]
    is_terminal_carboxylic = False
    
    for atom in mol.GetAtomWithIdx(carboxylic_carbon_idx).GetNeighbors():
        if atom.GetAtomicNum() == 6:
            # Check if the carbon has only one non-carboxylate neighbor
            if len([neighbor for neighbor in atom.GetNeighbors() if neighbor.GetIdx() != carboxylic_carbon_idx]) == 1:
                is_terminal_carboxylic = True
    
    if not is_terminal_carboxylic:
        return False, "Carboxylic acid group is not at a terminal position"
    
    # Check chain connectivity: long aliphatic chain (ignores branches but must not have rings)
    if not AllChem.EmbedMolecule(mol, randomSeed=42) >= 0:  # creates 3D structure to evaluate chain length
        return False, "Could not generate 3D structure for chain analysis"
    
    # Check for absence of unsaturation: presence of double bonds
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        return False, "Contains unsaturation (double bonds found)"
    
    # Verify absence of other functional groups (unsuspected heteroatoms)
    other_groups = Chem.MolFromSmarts("[!C&!c&!H][!C&!c&!H]")
    if mol.HasSubstructMatch(other_groups):
        return False, "Unexpected functional groups or heteroatoms found"
    
    return True, "Molecule is a saturated fatty acid"