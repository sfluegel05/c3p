"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is characterized by a deprotonated carboxylic acid group 
    and a long hydrocarbon chain, possibly with some functionalization.

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
    
    # Look for carboxylate group pattern [CX3](=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found (required for fatty acid anion)"
    
    # Check for a long hydrocarbon chain
    def has_sufficient_hydrocarbon_content(mol):
        carbo_chains = [len(match) for match in Chem.FindAllPathsOfLengthN(mol, 12, useBonds=False, fromAtoms=[atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])]
        return max(carbo_chains, default=0) >= 10  # At least 10 carbons considered as part of a longer chain

    if not has_sufficient_hydrocarbon_content(mol):
        return False, "Longest carbon chain is too short to be a fatty acid anion"

    # The overall flexibility to nuance functionalization over raw chain length
    return True, "Contains a carboxylate group and enough hydrocarbon content characteristic of a fatty acid anion"