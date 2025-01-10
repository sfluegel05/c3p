"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Classifies a molecule as a very long-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Convert the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Substructure pattern for CoA (coenzyme A) part
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"

    # Find the longest carbon chain
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    longest_chain_length = 0
    for atom in carbon_atoms:
        # Using a breadth-first search to find the longest path from each carbon atom
        queue = [(atom, 0)]
        visited = set([atom.GetIdx()])
        while queue:
            current_atom, length = queue.pop(0)
            if length > longest_chain_length:
                longest_chain_length = length
            
            for neighbor in current_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    queue.append((neighbor, length + 1))
                    visited.add(neighbor.GetIdx())
    
    if longest_chain_length <= 22:
        return False, f"Longest carbon chain is {longest_chain_length}, must be greater than 22"

    return True, "Contains CoA structure with a fatty acyl chain longer than C22"