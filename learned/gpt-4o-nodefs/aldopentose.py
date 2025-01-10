"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a monosaccharide with five carbon atoms and either an aldehyde functional group
    in its open-chain form or exists as a furanose or pyranose ring form.

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

    # Check for aldehyde group in open-chain form
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH2][CH](O)[CH](O)[CH2](O)")
    
    # Check for pyranose form (6-membered ring with oxygen, typical in sugars)
    pyranose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]1")
    
    # Check for furanose form (5-membered ring with oxygen, typical in sugars)
    furanose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)C1")
    
    has_structure = mol.HasSubstructMatch(aldehyde_pattern) or \
                    mol.HasSubstructMatch(pyranose_pattern) or \
                    mol.HasSubstructMatch(furanose_pattern)
    
    if not has_structure:
        return False, "No suitable aldehyde, pyranose, or furanose structure detected"
    
    # Check for sufficient oxygen atoms needed for hydroxyl groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:  # Expect at least 4 oxygens: 3 in OH groups + 1 in C=O/C-O
        return False, f"Insufficient oxygen atoms: expected at least 4, found {o_count}"
    
    # Exclude molecules with undesired functional groups
    undesired_patterns = [
        Chem.MolFromSmarts("P(=O)(O)O"),  # Phosphate esters
        Chem.MolFromSmarts("[CX3](=O)[O][H]")  # Carboxylic acids
    ]
    for pattern in undesired_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains undesired functional group"

    return True, "Molecule matches the characteristics of an aldopentose"