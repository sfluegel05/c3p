"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: CHEBI:29017 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative where one or more oxygens or hydroxyl groups
    are replaced by sulfur or -SR, where R can be hydrogen or any group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains any sulfur atoms
    has_sulfur = any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())
    if not has_sulfur:
        return False, "Molecule does not contain sulfur"

    # Check if the molecule contains a carbohydrate backbone
    carbohydrate_pattern = Chem.MolFromSmarts("[OX2]([CX4])[CX4][CX4]")
    carbohydrate_matches = mol.GetSubstructMatches(carbohydrate_pattern)
    if not carbohydrate_matches:
        return False, "No carbohydrate backbone found"

    # Check if any sulfur atom is attached to the carbohydrate backbone
    for sulfur_idx in [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]:
        sulfur = mol.GetAtomWithIdx(sulfur_idx)
        neighbors = [mol.GetAtomWithIdx(nbr_idx) for nbr_idx in sulfur.GetNeighbors()]
        if any(nbr.GetAtomicNum() == 8 and nbr.IsInRing() for nbr in neighbors):
            return True, "Sulfur atom replacing oxygen or hydroxyl group in carbohydrate backbone"

    return False, "Sulfur atom not replacing oxygen or hydroxyl group in carbohydrate backbone"