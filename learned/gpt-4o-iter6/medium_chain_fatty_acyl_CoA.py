"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for thioester linkage (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage to CoA found"
    
    # Full Coenzyme A moiety pattern including thioester linkage and enough neighboring atoms
    full_coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(full_coa_pattern):
        return False, "Full CoA structure not correctly matched"
    
    # Calculate and check the length of the carbon chain bound to the thioester group
    fatty_chain_length = get_fatty_acid_chain_length(mol)
    if not (6 <= fatty_chain_length <= 12):
        return False, f"Aliphatic chain length of {fatty_chain_length} not within medium-chain range (6-12 carbons)"
    
    return True, "Molecule is a medium-chain fatty acyl-CoA with proper CoA moiety and chain length"

def get_fatty_acid_chain_length(mol):
    """
    Utility to count the longest continuous carbon chain starting from the thioester linkage.
    """
    chain_lengths = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() >= 2:
            # Look for paths starting at the carbonyl carbon (COC in COS)
            paths = Chem.rdmolops.GetShortestPaths(mol, atom.GetIdx())
            chain_lengths.extend(max(len(path) for path in paths if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 1 for idx in path)))
    return max(chain_lengths) if chain_lengths else 0

# Example test
example_smiles = "CCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
print(is_medium_chain_fatty_acyl_CoA(example_smiles))