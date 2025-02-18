"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid (CHEBI: ???)
A methyl-branched fatty acid is a branched-chain fatty acid containing methyl branches only.
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid has a carboxylic acid group and only methyl branches.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid = MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Get indices of carboxylic acid carbon and oxygen atoms
    cooh_matches = mol.GetSubstructMatches(carboxylic_acid)
    cooh_carbons = {match[0] for match in cooh_matches}  # Carbonyl carbon indices

    # Check all non-COOH carbons for non-methyl branches
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() not in cooh_carbons:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in cooh_carbons:
                    # Check if neighbor has other carbon connections besides current atom
                    neighbor_carbon_count = 0
                    for n in neighbor.GetNeighbors():
                        if n.GetAtomicNum() == 6 and n != atom:
                            neighbor_carbon_count += 1
                    if neighbor_carbon_count > 0:
                        return False, f"Non-methyl branch at carbon {neighbor.GetIdx()}"

    # Verify there's at least one methyl branch (to exclude linear fatty acids)
    methyl_pattern = MolFromSmarts("[CH3]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    has_methyl_branch = False
    for match in methyl_matches:
        methyl_carbon = mol.GetAtomWithIdx(match[0])
        # Check if methyl is attached to a non-COOH carbon
        for neighbor in methyl_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in cooh_carbons:
                has_methyl_branch = True
                break
        if has_methyl_branch:
            break
    if not has_methyl_branch:
        return False, "No methyl branches found"

    return True, "Methyl-branched fatty acid with only methyl substituents"