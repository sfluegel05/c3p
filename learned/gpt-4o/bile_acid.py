"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is a hydroxy-5beta-cholanic acid occurring in bile.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update substructure pattern for cholanic acid backbone with 5beta configuration
    cholanic_base_pattern = Chem.MolFromSmarts('[CH]1([C@@]2([C@@H](CC[C@]3([C@@H](CCC[C@@H]4[C@H]3CC[C@@H]4O)C2)C)C1)O)')

    if not mol.HasSubstructMatch(cholanic_base_pattern):
        return False, "No cholanic acid backbone with 5β configuration"

    # Count the number of hydroxy groups
    hydroxy_group = Chem.MolFromSmarts('[CX4][OH]')
    ho_count = len(mol.GetSubstructMatches(hydroxy_group))

    if ho_count < 1:
        return False, "No hydroxy groups identified"

    # Presence of carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    return True, "Contains cholanic acid backbone with 5β configuration and typical functional groups of bile acids"