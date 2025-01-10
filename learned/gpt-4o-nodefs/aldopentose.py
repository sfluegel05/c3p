"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a monosaccharide with five carbon atoms and either an aldehyde functional group
    in its open-chain form or exists as a cyclic form like furanose or pyranose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly 5 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected 5 carbon atoms, found {c_count}"

    # Check for aldehyde group: pattern C(=O) connected to a terminal H (open-chain form)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)CO")
    # Also check for cyclic ring forms (furanose/pyranose)
    furanose_pyranose_pattern = Chem.MolFromSmarts("C1OC(O)C(CO)O1")
    has_aldopentose_structure = mol.HasSubstructMatch(aldehyde_pattern) or mol.HasSubstructMatch(furanose_pyranose_pattern)
    
    if not has_aldopentose_structure:
        return False, "No aldehyde group or furanose/pyranose ring detected"
    
    # Check for sufficient oxygen atoms needed for hydroxyl groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Insufficient oxygen atoms: expected at least 4, found {o_count}"
    
    # Check for undesired groups like phosphates or carboxylic acids
    undesired_group_patterns = [
        Chem.MolFromSmarts("P(=O)(O)O"),
        Chem.MolFromSmarts("C(=O)[OH]"),
    ]
    for pattern in undesired_group_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, f"Contains undesired functional group"

    return True, "Molecule matches the characteristics of an aldopentose"