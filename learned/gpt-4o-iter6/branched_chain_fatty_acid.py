"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid has one or more alkyl substituents and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (R-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Detailed carbon counting for chain identification
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:  # Relax threshold, typical for smaller branched acids
        return False, "Insufficient carbon count for alkyl chain"

    # Look for branching - general multiple substituent model
    branch_pattern = Chem.MolFromSmarts("[CX4][CX4,CX3]([CX4,H])[CX4,CX3,H]")  # Capture different branching forms
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No branch points found"

    # Filter out some complex peptide-like or too aromatic structures
    # This is a simple check and can be refined more in real scenarios
    peptide_like_multiplier = mol.GetRingInfo().NumRings() + mol.GetNumAtoms() - c_count
    if peptide_like_multiplier > 10:
        return False, "Too complex or peptide-like structure for a simple alkyl branch"

    return True, "Contains a carbon chain with branch points and a carboxylic acid group with acceptable complexity"