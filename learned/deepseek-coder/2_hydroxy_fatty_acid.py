"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: 2-hydroxy fatty acid
Definition: Any fatty acid with a hydroxy functional group in the alpha- or 2-position.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a carboxylic acid group and a hydroxyl group at the 2-position,
    along with a sufficiently long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Use a more flexible SMARTS pattern to detect the hydroxyl group at the 2-position
    # relative to the carboxylic acid group
    # This pattern allows for branching and double bonds in the carbon chain
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1][CX4H1,CX3H1][OX2H1]")
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No hydroxyl group at the 2-position relative to the carboxylic acid group"

    # Check for a sufficiently long carbon chain (fatty acid)
    # Allow for shorter chains, as some 2-hydroxy fatty acids can have fewer than 6 carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, "Carbon chain too short to be a fatty acid"

    # Additional check to ensure the molecule is a fatty acid (linear or branched carbon chain)
    # Exclude molecules with complex structures like rings or multiple functional groups
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Molecule contains rings, not a typical fatty acid"

    return True, "Contains a carboxylic acid group, a hydroxyl group at the 2-position, and a sufficiently long carbon chain"