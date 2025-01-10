"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is characterized by one or more amino groups replacing hydroxyl groups on a sugar molecule.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for amino group which should replace a hydroxyl group
    amino_group_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]")  # Any amine (primary, secondary, tertiary)

    # Check we have an amino group
    if not mol.HasSubstructMatch(amino_group_pattern):
        return False, "No amino groups found"

    # Basic pattern for a sugar-like structure, matches pyranose/furanose rings commonly
    pyranose_pattern = Chem.MolFromSmarts("C1OC([C@H](O)[C@@H](O)[C@H](O))C(O)C1")
    furanose_pattern = Chem.MolFromSmarts("C1OC([C@H](O)[C@@H](O)[C@@H](O)C1)O")
    
    # Check for sugar-like structure
    has_sugar_ring = mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern)
    if not has_sugar_ring:
        return False, "No sugar-like structure found"

    # Check replacements of the hydroxyl groups by amine groups
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            # Inspect bonded atoms to check if it's possible to replace OH
            nitrogen_neighbors = [a.GetAtomicNum() for a in atom.GetNeighbors()]
            if 6 in nitrogen_neighbors:  # Checking for carbon around the nitrogen
                neighboring_carbon = next(a for a in atom.GetNeighbors() if a.GetSymbol() == 'C')
                # If a carbon atom's hydroxyl is replaced by amine so check removal of oxygen which should be adjacent
                carbon_oxygen_bond = [a.GetAtomicNum() for a in neighboring_carbon.GetNeighbors()]
                if 8 not in carbon_oxygen_bond:
                    return True, "Contains a sugar backbone with amino groups substituting hydroxyl groups."

    return False, "No correct hydroxyl replacement was found by an amino group."