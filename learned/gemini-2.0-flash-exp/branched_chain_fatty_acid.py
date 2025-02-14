"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid has a carboxylic acid group, a long carbon chain and one or more alkyl substituents.

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

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for at least 5 carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, "Too few carbons for a fatty acid"

    # Check for branching point (carbon connected to at least 3 other non-hydrogen atoms)
    branch_pattern = Chem.MolFromSmarts("C([CX4])([CX4])[CX4]")
    if not mol.HasSubstructMatch(branch_pattern):
        branch_pattern = Chem.MolFromSmarts("C([CX4])([CX4])([CX3,CX2])") #allows for double/triple bonds
        if not mol.HasSubstructMatch(branch_pattern):
            return False, "No branching point found"

    #check for long chain - add minimum number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Chains too short to be a fatty acid"
    
    return True, "Contains a carboxylic acid group, a long chain and one or more branching points"