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

    # Check for branching point (carbon connected to at least 3 other non-hydrogen atoms)
    branch_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6: # Check only carbon atoms
             degree = len([neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() != 1]) #exclude H
             if degree >= 3:
                branch_found = True
                break

    if not branch_found:
        return False, "No branching point found"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 130: # Fatty acids usually have > 6 carbons (minimum MW ~130)
         return False, "Molecular weight too low for fatty acid"
    
    return True, "Contains a carboxylic acid group, a long chain and one or more branching points"