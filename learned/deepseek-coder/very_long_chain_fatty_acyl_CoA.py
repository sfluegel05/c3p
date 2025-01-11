"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:XXXXX very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA is defined as a fatty acyl-CoA with a fatty acyl chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA moiety pattern (simplified)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Find the fatty acyl chain
    # Look for the carbonyl group attached to the CoA sulfur
    fatty_acid_start = Chem.MolFromSmarts("[CX3](=O)SCCNC(=O)CCNC(=O)")
    matches = mol.GetSubstructMatches(fatty_acid_start)
    if not matches:
        return False, "No fatty acyl chain found attached to CoA"

    # Get the atom connected to the carbonyl carbon
    carbonyl_carbon = matches[0][0]
    fatty_chain_atom = mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors()[0].GetIdx()

    # Traverse the chain and count carbons
    chain_length = 0
    visited = set()
    stack = [fatty_chain_atom]
    
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        atom = mol.GetAtomWithIdx(current)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_length += 1
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() in [1, 6]:  # Hydrogen or carbon
                    stack.append(neighbor.GetIdx())

    # Subtract 1 for the carbonyl carbon
    chain_length -= 1

    if chain_length > 22:
        return True, f"Fatty acyl chain length is {chain_length} (greater than C22)"
    else:
        return False, f"Fatty acyl chain length is {chain_length} (not greater than C22)"