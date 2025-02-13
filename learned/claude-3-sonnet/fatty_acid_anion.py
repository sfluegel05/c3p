"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: CHEBI:36195 fatty acid anion
The conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group of the corresponding fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylate group (-COO-)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Look for carbon chain of variable length (at least 3 carbons)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~*")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No carbon chain found"

    # Exclude specific substructures not found in typical fatty acid anions
    exclude_patterns = [
        Chem.MolFromSmarts("[#7r5,#7r6,#7r7,#7r8,#7r9]"),  # Exclude complex heterocycles
        Chem.MolFromSmarts("[#16]"),  # Exclude sulfur atoms
        Chem.MolFromSmarts("[#7,#8]~[#6]~[#7,#8]"),  # Exclude amide groups
    ]
    for pattern in exclude_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains substructures not typical of fatty acid anions"

    # Check molecular weight - fatty acid anions typically > 150 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for fatty acid anion"

    # Optionally, handle stereochemistry if needed
    # ...

    return True, "Contains a carboxylate group and a carbon chain, consistent with a fatty acid anion"