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
    
    # Check for thioester pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage to CoA found"
    
    # Broad Coenzyme A moiety pattern including adenine-related parts
    coa_moiety_pattern = Chem.MolFromSmarts("NCC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(O)OC[C@H]1O[C@H](O)[C@@H]1OP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_moiety_pattern):
        return False, "CoA moiety structure not matched"
    
    # Calculate the fatty acid chain length bound to the thioester group
    fatty_chain_length = calculate_chain_length(mol, thioester_pattern)

    # Medium-chain fatty acid range is about 6-12 carbons
    if not (6 <= fatty_chain_length <= 12):
        return False, f"Aliphatic chain length of {fatty_chain_length} not within medium-chain range (6-12 carbons)"
    
    return True, "Molecule is a medium-chain fatty acyl-CoA with proper CoA moiety and chain length"


def calculate_chain_length(mol, pattern):
    """
    Uses depth-first search to calculate the longest carbon chain length starting from the thioester carbon.

    Args:
        mol: molecule object obtained from RDKit.
        pattern: SMARTS pattern that identifies thioester linkage.

    Returns:
        Longest chain length found starting at thioester carbon.
    """
    # Find the carbon linked to sulfur in thioester
    match_idx = mol.GetSubstructMatch(pattern)
    if not match_idx:
        return 0  # Failed to find starting point

    start_atom = match_idx[1]  # S carbon in C(=O)SCC

    # Perform a DFS starting with the carbon in "-COS"
    def dfs(atom_idx, seen_atoms):
        seen_atoms.add(atom_idx)
        max_chain_length = 1

        for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
            n_idx = neighbor.GetIdx()
            # Avoid revisiting atoms, and only traverse carbon atoms
            if n_idx not in seen_atoms and neighbor.GetAtomicNum() == 6:
                max_chain_length = max(max_chain_length, 1 + dfs(n_idx, seen_atoms))

        seen_atoms.remove(atom_idx)
        return max_chain_length

    return dfs(start_atom, set())

# Example test
example_smiles = "CCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
print(is_medium_chain_fatty_acyl_CoA(example_smiles))