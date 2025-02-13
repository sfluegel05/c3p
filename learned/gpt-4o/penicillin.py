"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

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
    # This should encompass a thiazolidine ring fused to a beta-lactam
    penam_core_pattern = Chem.MolFromSmarts("C1C(N2C(S1)C(C2=O)N)C")
    if not mol.HasSubstructMatch(penam_core_pattern):
        return False, "Penam core structure not found"

    # Check for two methyl groups at position 2 of thiazolidine ring
    methyl_groups = Chem.MolFromSmarts("C1(C)C(S1)(C)")
    if not mol.HasSubstructMatch(methyl_groups):
        return False, "Two methyl groups not found at position 2"

    # Check for carboxylate at position 3
    carboxylate_group = Chem.MolFromSmarts("C(=O)[O-]")
    carboxyl_group = Chem.MolFromSmarts("C(=O)O")  # For carboxylic acid form
    if not mol.HasSubstructMatch(carboxylate_group) and not mol.HasSubstructMatch(carboxyl_group):
        return False, "Carboxylate group not found at position 3"

    # Check for carboxamido group at position 6
    carboxamido_group = Chem.MolFromSmarts("C(=O)NC")
    if not mol.HasSubstructMatch(carboxamido_group):
        return False, "Carboxamido group not found at position 6"

    # Ensure stereochemistry matches penicillin, focusing on key atoms
    is_stereochemical_correct = False
    for atom in mol.GetAtoms():
        if atom.GetChiralTag() not in [Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW]:
            continue
        # Rough check on stereocenter: atom must be part of penam structure
        if any(neigh.IsInRing() and neigh.GetSymbol() == 'S' for neigh in atom.GetNeighbors()):
            is_stereochemical_correct = True
            break

    if not is_stereochemical_correct:
        return False, "Stereochemistry does not match penicillin"

    return True, "Molecule matches the penicillin structure"

# Example test
example_smiles = '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)CCC[C@H](N)C(O)=O)C(O)=O'
print(is_penicillin(example_smiles))