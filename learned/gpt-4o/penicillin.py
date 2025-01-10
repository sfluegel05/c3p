"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin must have a specific penam core structure with additional specific substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define penam core SMARTS pattern
    penam_core_pattern = Chem.MolFromSmarts("C1[C@H]2S[C@](C)(C)N2C(=O)[C@H]1NC(=O)")
    if not mol.HasSubstructMatch(penam_core_pattern):
        return False, "Penam core structure not found"

    # Check for two methyl groups at position 2 of thiazolidine ring
    methyl_group_pattern = Chem.MolFromSmarts("C1(C)(C)S[C@@H]1C")
    if not mol.HasSubstructMatch(methyl_group_pattern):
        return False, "Two methyl groups not found at position 2"

    # Check for carboxylate at position 3
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    carboxyl_group_pattern = Chem.MolFromSmarts("C(=O)O")  # Carboxylic acid form
    if not mol.HasSubstructMatch(carboxylate_pattern) and not mol.HasSubstructMatch(carboxyl_group_pattern):
        return False, "Carboxylate group not found at position 3"

    # Check for carboxamido group at position 6
    carboxamido_pattern = Chem.MolFromSmarts("C(=O)NC")
    if not mol.HasSubstructMatch(carboxamido_pattern):
        return False, "Carboxamido group not found at position 6"

    # Check stereochemistry matches typical penicillin
    stereochemistry_correct = False
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() not in [Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW]:
            continue
        # Check that this atom is in a ring containing S, matching expected penicillin configurations
        if any(neigh.IsInRing() and neigh.GetSymbol() == 'S' for neigh in atom.GetNeighbors()):
            stereochemistry_correct = True
            break

    if not stereochemistry_correct:
        return False, "Stereochemistry does not match penicillin"

    return True, "Molecule matches the penicillin structure"