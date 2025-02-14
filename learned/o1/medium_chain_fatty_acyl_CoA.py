"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:60903 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA results from the condensation of coenzyme A with a medium-chain fatty acid.
    Medium-chain fatty acids are aliphatic chains with 6 to 13 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define coenzyme A fragment SMARTS
    coa_smarts = '[C@H](OP(=O)(O)OCC1OC(C(O)C1OP(=O)(O)O)n2cnc3c(N)ncnc23)O'  # Simplified CoA part

    # Define fatty acyl-CoA thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts('C(=O)SCCNC(=O)CCNC(=O)' + coa_smarts)
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No fatty acyl-CoA structure detected"

    # Find the thioester carbonyl carbon
    thioester_pattern = Chem.MolFromSmarts('C(=O)SC')
    matches = mol.GetSubstructMatches(thioester_pattern)
    if not matches:
        return False, "Thioester linkage not found"

    # Get the carbonyl carbon atom index in the thioester
    carbonyl_c_index = matches[0][0]

    # Traverse the acyl chain from the carbonyl carbon
    carbons_in_chain = set()
    queue = [mol.GetAtomWithIdx(carbonyl_c_index)]
    visited = set()
    contains_aromatic = False

    while queue:
        atom = queue.pop()
        atom_idx = atom.GetIdx()

        if atom_idx in visited:
            continue
        visited.add(atom_idx)

        if atom.GetAtomicNum() == 6:
            carbons_in_chain.add(atom_idx)
            if atom.GetIsAromatic():
                contains_aromatic = True

        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in visited:
                continue
            bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
            # Avoid crossing the thioester bond towards CoA
            if bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetAtomicNum() != 16:
                queue.append(neighbor)

    # Exclude aromatic acyl chains
    if contains_aromatic:
        return False, "Aromatic acyl chains are not considered fatty acids"

    carbon_count = len(carbons_in_chain)

    # Medium-chain fatty acids have 6 to 13 carbons
    if 6 <= carbon_count <= 13:
        return True, f"Contains fatty acyl chain with {carbon_count} carbons (medium-chain)"
    else:
        return False, f"Fatty acyl chain with {carbon_count} carbons is not medium-chain"