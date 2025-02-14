"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: CHEBI:52634 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine is a ceramide consisting of sphinganine where one of the amino hydrogens is substituted by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define amide bond pattern (N-C(=O)-C)
    amide_pattern = Chem.MolFromSmarts('N-C(=O)-C')
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond (N-acyl group) found"

    # Assume the nitrogen in the amide is the sphinganine nitrogen
    # Check for sphinganine backbone connected to this nitrogen
    for match in amide_matches:
        amide_nitrogen_idx = match[0]
        amide_nitrogen = mol.GetAtomWithIdx(amide_nitrogen_idx)

        # Get the alpha carbon (connected to the nitrogen but not the carbonyl carbon)
        neighbors = [nbr for nbr in amide_nitrogen.GetNeighbors() if nbr.GetIdx() != match[1]]
        if not neighbors:
            continue
        alpha_carbon = neighbors[0]

        # Check if alpha carbon has a hydroxyl group
        has_alpha_hydroxyl = False
        for nbr in alpha_carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # Oxygen
                has_alpha_hydroxyl = True
                break
        if not has_alpha_hydroxyl:
            continue

        # Get beta carbon (next carbon in the chain)
        beta_carbons = [nbr for nbr in alpha_carbon.GetNeighbors() if nbr.GetIdx() != amide_nitrogen.GetIdx()]
        if not beta_carbons:
            continue
        beta_carbon = beta_carbons[0]

        # Get gamma carbon
        gamma_carbons = [nbr for nbr in beta_carbon.GetNeighbors() if nbr.GetIdx() != alpha_carbon.GetIdx()]
        if not gamma_carbons:
            continue
        gamma_carbon = gamma_carbons[0]

        # Check if gamma carbon has a hydroxyl group
        has_gamma_hydroxyl = False
        for nbr in gamma_carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # Oxygen
                has_gamma_hydroxyl = True
                break
        if not has_gamma_hydroxyl:
            continue

        # Check for long aliphatic chain from gamma carbon
        chain_length = 0
        visited = set()
        to_visit = [gamma_carbon]
        while to_visit:
            current_atom = to_visit.pop()
            if current_atom.GetIdx() in visited:
                continue
            visited.add(current_atom.GetIdx())
            if current_atom.GetAtomicNum() == 6 and current_atom.GetDegree() <= 4:
                chain_length += 1
                for nbr in current_atom.GetNeighbors():
                    if nbr.GetIdx() not in visited and nbr.GetAtomicNum() == 6:
                        to_visit.append(nbr)

        if chain_length < 12:
            continue  # Chain is too short to be sphinganine

        # Check acyl chain length (from carbonyl carbon)
        carbonyl_carbon = mol.GetAtomWithIdx(match[1])
        acyl_chain_length = 0
        visited = set()
        to_visit = [nbr for nbr in carbonyl_carbon.GetNeighbors() if nbr.GetIdx() != amide_nitrogen_idx]
        while to_visit:
            current_atom = to_visit.pop()
            if current_atom.GetIdx() in visited:
                continue
            visited.add(current_atom.GetIdx())
            if current_atom.GetAtomicNum() == 6 and current_atom.GetDegree() <= 4:
                acyl_chain_length += 1
                for nbr in current_atom.GetNeighbors():
                    if nbr.GetIdx() not in visited and nbr.GetAtomicNum() == 6:
                        to_visit.append(nbr)

        if acyl_chain_length < 12:
            continue  # Acyl chain is too short

        return True, "Contains sphinganine backbone with N-acyl group forming an amide bond"

    return False, "No sphinganine backbone with N-acyl group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:52634',
        'name': 'N-acylsphinganine',
        'definition': "A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.",
        'parents': ['CHEBI:33598', 'CHEBI:33709']
    },
    'config': {}
}