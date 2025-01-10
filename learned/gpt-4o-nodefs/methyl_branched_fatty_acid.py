"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid or carboxylate group -COOH or -COO- (smarts pattern: C(=O)[OH,[O-]])
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O,OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid or carboxylate group found"

    # Look for at least one methyl group attached to a secondary or tertiary carbon
    methyl_branch_pattern = Chem.MolFromSmarts("[CX4](C)([CH3])[CX4]")
    if not mol.HasSubstructMatch(methyl_branch_pattern):
        return False, "No methyl branching found"

    # Count aliphatic carbons to ensure sufficient chain length for fatty acids
    aliphatic_carbon_count = sum(1 for atom in mol.GetAtoms() 
                                 if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic())
    if aliphatic_carbon_count < 8:
        return False, "Not enough carbon atoms for a fatty acid backbone"

    return True, "Contains carboxylic acid or carboxylate group and methyl branching in chain"