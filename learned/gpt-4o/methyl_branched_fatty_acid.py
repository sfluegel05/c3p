"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.

    A methyl-branched fatty acid is defined as any branched-chain fatty acid
    containing methyl branches only and having a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check if every branch is a methyl group
    # A methyl branch is a terminal carbon bonded to another carbon (no ring), attached to the backbone
    methyl_branch_pattern = Chem.MolFromSmarts("[C;X4;H3]")
    methyl_branches = mol.GetSubstructMatches(methyl_branch_pattern)
    if not methyl_branches:
        return False, "No methyl branches found"

    # Ensure branches are on the fatty acid backbone (should not form multiple branching beyond methyl)
    for atom in mol.GetAtoms():
        # Check if the central carbon has any non-methyl groups attached erroneously
        if atom.GetSymbol() == 'C' and sum(1 for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'C') > 1:
            nbr_symbols = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            # More than one neighbor carbon usually indicates a branch (must ensure neighbors are methyl only)
            if sum(symbol == 'C' for symbol in nbr_symbols) - 1 > sum(atom.GetSymbol() == 'H' and atom.GetDegree() == 1 for atom in atom.GetNeighbors()):
                return False, "Non-methyl branching found"
    
    # Estimate chain length must encompass carbon counts typical of fatty acids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5 or c_count > 40:  # Extended to a limit based on known examples
        return False, f"Carbon chain length off typical range for fatty acids: {c_count} carbons"

    return True, "Molecule is a methyl-branched fatty acid"

# Usage Example
# result, reason = is_methyl_branched_fatty_acid("OC(=O)C(CCCC(C)C)C")
# print(result, reason)