"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA typically features an acyl chain with at least 22 carbons
    attached to a CoA moiety via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the thioester group pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester group found"

    # Identify the CoA moiety using a more complete SMARTS pattern, including adenylate moiety
    coa_pattern = Chem.MolFromSmarts("[P](=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not detected"

    # Count carbon atoms in the longest aliphatic chain prior to the thioester linkage
    # Assume the aliphatic chain ends at the carbonyl carbon in thioester
    fragments = mol.GetSubstructMatches(thioester_pattern)
    longest_chain_length = 0
    if fragments:
        for frag in fragments:
            carbon_chain = Chem.RWMol(mol).GetFragmentAtoms(frag[0])
            carbon_count = sum(1 for atom in carbon_chain if atom.GetAtomicNum() == 6)
            longest_chain_length = max(longest_chain_length, carbon_count)
    
    if longest_chain_length < 22:
        return False, "Aliphatic chain not sufficiently long (very long-chain fatty acids typically have ≥22 carbons)"

    return True, "Contains a long fatty acyl chain and CoA moiety with thioester bond"

# Example test
example_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(=O)O)n1cnc2c(N)ncnc12"
result, reason = is_very_long_chain_fatty_acyl_CoA(example_smiles)
print(result, reason)