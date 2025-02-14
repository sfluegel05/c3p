"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:57373 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA results from the formal condensation of the thiol group
    of coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define CoA SMARTS pattern (simplified for matching)
    coa_smarts = """
    NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](COP(O)(O)=O)[C@H](O)[C@@H]1OP(O)(O)=O
    """
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Invalid CoA SMARTS pattern"

    # Find CoA moiety
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "No CoA moiety found"

    # Define thioester pattern C(=O)S
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"

    # Find the thioester carbon and sulfur atoms
    for match in thioester_matches:
        c_idx = match[0]  # Carbonyl carbon index
        s_idx = match[2]  # Sulfur index

        # Check if sulfur is connected to CoA
        is_sulfur_in_coa = False
        for coa_match in coa_matches:
            if s_idx in coa_match:
                is_sulfur_in_coa = True
                break
        if not is_sulfur_in_coa:
            continue  # Thioester not connected to CoA

        # Identify acyl chain starting from carbonyl carbon
        carbonyl_carbon = mol.GetAtomWithIdx(c_idx)

        # Exclude bonds to oxygen and sulfur
        acyl_neighbors = [
            nbr for nbr in carbonyl_carbon.GetNeighbors()
            if nbr.GetAtomicNum() != 8 and nbr.GetAtomicNum() != 16
        ]
        if not acyl_neighbors:
            continue  # No acyl chain attached to carbonyl carbon
        acyl_start_atom = acyl_neighbors[0]

        # Use BFS traversal to get acyl chain atoms
        acyl_chain_atoms = set()
        to_visit = [acyl_start_atom.GetIdx()]
        while to_visit:
            current_idx = to_visit.pop()
            current_atom = mol.GetAtomWithIdx(current_idx)
            if current_atom.GetAtomicNum() in (1, 8, 15, 16, 7):  # Exclude non-carbon atoms
                continue
            if current_idx in acyl_chain_atoms:
                continue
            acyl_chain_atoms.add(current_idx)
            # Add neighbor carbons to visit
            for neighbor in current_atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if neighbor.GetAtomicNum() == 6 and nbr_idx != carbonyl_carbon.GetIdx():
                    to_visit.append(nbr_idx)

        # Count the number of carbons in the acyl chain
        carbon_count = len(acyl_chain_atoms) + 1  # Include carbonyl carbon

        if 13 <= carbon_count <= 22:
            return True, f"Contains CoA moiety linked via thioester to a fatty acyl chain of length {carbon_count} carbons"
        else:
            return False, f"Acyl chain length is {carbon_count}, not in the range 13-22 carbons"

    return False, "Thioester group not connected to CoA moiety or acyl chain length not in range"