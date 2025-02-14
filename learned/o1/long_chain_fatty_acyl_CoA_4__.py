"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:77989 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    Checks for the presence of a CoA moiety, a thioester linkage, and a long acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple:
            bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
            str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize molecule to standardize bonding patterns
    try:
        Chem.Kekulize(mol)
    except Chem.KekulizeException:
        pass  # Ignore kekulization errors

    # Define a comprehensive CoA SMARTS pattern
    coa_smarts = """
    NC(=O)CCNC(=O)[C@@H](O)[C@](C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H](COP(=O)([O-])[O-])[C@@H](O)[C@H]1O
    n2cnc3c(N)ncnc23
    """
    coa_smarts = coa_smarts.replace('\n', '').replace(' ', '')
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Invalid CoA SMARTS pattern"

    # Check for CoA moiety
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define thioester linkage pattern: C(=O)S
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Find the carbonyl carbon atom of the thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"
    # Assume the first match is the relevant one
    thioester_carb_idx = thioester_matches[0][0]
    thioester_sulfur_idx = thioester_matches[0][2]

    # Set of atom indices belonging to CoA moiety
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    coa_atom_indices = set()
    for match in coa_matches:
        coa_atom_indices.update(match)

    # Trace the acyl chain starting from the carbonyl carbon
    acyl_chain_atoms = set()
    atoms_to_visit = [thioester_carb_idx]
    while atoms_to_visit:
        current_atom_idx = atoms_to_visit.pop()
        if current_atom_idx in acyl_chain_atoms:
            continue
        if current_atom_idx in coa_atom_indices:
            continue
        acyl_chain_atoms.add(current_atom_idx)
        current_atom = mol.GetAtomWithIdx(current_atom_idx)
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in acyl_chain_atoms:
                continue
            # Stop at sulfur atom to prevent crossing into CoA moiety
            if neighbor_idx == thioester_sulfur_idx:
                continue
            # Exclude atoms that are part of CoA
            if neighbor_idx in coa_atom_indices:
                continue
            atoms_to_visit.append(neighbor_idx)

    # Count the number of carbon atoms in the acyl chain, excluding the carbonyl carbon
    acyl_chain_carbons = [
        idx for idx in acyl_chain_atoms
        if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and idx != thioester_carb_idx
    ]
    acyl_chain_length = len(acyl_chain_carbons)

    if acyl_chain_length < 12:
        return False, f"Acyl chain too short ({acyl_chain_length} carbons), need at least 12"

    return True, f"Contains CoA moiety with long acyl chain of {acyl_chain_length} carbons"