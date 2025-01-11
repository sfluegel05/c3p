"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any medium-chain fatty acid.
"""

from rdkit import Chem

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

    # Identify the thioester linkage: C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Expected one thioester linkage, found {len(thioester_matches)}"

    # Get indices of carbonyl carbon and sulfur atom in the thioester linkage
    carbonyl_c_idx = thioester_matches[0][0]
    sulfur_idx = thioester_matches[0][2]

    # Collect atoms of the acyl chain excluding the sulfur atom
    visited_atoms = set()

    def traverse_acyl_chain(atom_idx, exclude_idx):
        """
        Recursively traverse the acyl chain starting from the carbonyl carbon,
        excluding the path towards the sulfur atom (CoA moiety).
        """
        if atom_idx in visited_atoms:
            return
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx != exclude_idx:
                traverse_acyl_chain(neighbor_idx, exclude_idx)

    # Start traversal from the carbonyl carbon, exclude sulfur atom
    traverse_acyl_chain(carbonyl_c_idx, sulfur_idx)

    # Count the number of carbons in the acyl chain
    num_carbons = sum(1 for idx in visited_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)

    # Check if the number of carbons is within medium-chain length (6-12)
    if num_carbons < 6:
        return False, f"Acyl chain too short ({num_carbons} carbons), not medium-chain"
    if num_carbons > 12:
        return False, f"Acyl chain too long ({num_carbons} carbons), not medium-chain"

    # Verify the presence of coenzyme A moiety
    # Define a SMARTS pattern for key features of CoA
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not found"

    return True, f"Contains medium-chain acyl group ({num_carbons} carbons) attached to coenzyme A via thioester linkage"