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
    
    # Look for terminal carboxylic acid or carboxylate group -COOH or -COO- (smarts pattern: [CH2,CH][CX3](=O)[OH])
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CH2,CH][CX3](=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid or carboxylate group found"

    # Look for methyl group branching not on terminal positions
    methyl_branch_pattern = Chem.MolFromSmarts("[CX4H2]C(C)(C)[C;!H0]")  # Methyl branch in chain
    if not mol.HasSubstructMatch(methyl_branch_pattern):
        return False, "No methyl branching found in the chain"

    # Count total carbon atoms to ensure sufficient chain length for fatty acids
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 10:  # Adjusted minimum length to ensure it's considered a fatty acid
        return False, "Not enough carbon atoms for a fatty acid backbone"

    return True, "Contains carboxylic acid or carboxylate group and appropriate methyl branching in chain"