"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA 
    based on its SMILES string. Very long-chain fatty acyl-CoAs 
    have fatty acyl groups with chain lengths greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA,
              False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify Coenzyme A part by locating substructures
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A (CoA) core structure found"
    
    # Locate the acyl chain part, i.e., part before thioester linkage
    thioester_idx = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]  # S atom in thioester linkage
    if not thioester_idx:
        return False, "No thioester linkage found"

    # Trace acyl chain: Start from any carbon directly attached to sulfur
    acyl_chain_atoms = set()
    visited = set()

    def atom_dfs(atom):
        if atom.GetIdx() in visited or atom.GetAtomicNum() != 6:  # Only consider carbons
            return
        acyl_chain_atoms.add(atom.GetIdx())
        visited.add(atom.GetIdx())
        # Check linear connections and avoid cyclic/branched paths
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in acyl_chain_atoms:
                atom_dfs(neighbor)

    thioester_sulfur = mol.GetAtomWithIdx(thioester_idx[0])
    for neighbor in thioester_sulfur.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Start DFS tree from carbon atoms
            atom_dfs(neighbor)
            break  # Consider only one primary chain extending from the attachment point

    # Count the number of carbon atoms in the acyl chain
    carbon_count = len(acyl_chain_atoms)
    
    # Determine if it's a very long-chain fatty acyl-CoA
    if carbon_count > 22:
        return True, f"Fatty acyl chain length is {carbon_count} (> 22 carbons, valid)"

    return False, f"Fatty acyl chain length is {carbon_count} (<= 22 carbons, not valid)"