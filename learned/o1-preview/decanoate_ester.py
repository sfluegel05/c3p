"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:36027 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is a fatty acid ester resulting from the formal condensation of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester group pattern
    ester_pattern = Chem.MolFromSmarts('C(=O)O[!$([O-])]')  # Matches ester functional group

    # Find ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check for decanoyl chain
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon atom index
        # Start traversal from carbonyl carbon
        visited = set()
        found_decanoyl = False

        def traverse(atom_idx, depth):
            nonlocal found_decanoyl
            if depth > 10:
                return
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                return  # Only traverse carbon atoms
            visited.add(atom_idx)
            if atom.GetDegree() == 1 and depth == 10:
                # Found a chain of length 10 carbons (including carbonyl carbon)
                found_decanoyl = True
                return
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
                if neighbor_idx not in visited and bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetAtomicNum() == 6:
                    traverse(neighbor_idx, depth+1)
            visited.remove(atom_idx)

        traverse(carbonyl_c_idx, 1)
        if found_decanoyl:
            return True, "Contains decanoate ester group"

    return False, "No decanoate ester group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36027',
        'name': 'decanoate ester',
        'definition': 'A fatty acid ester resulting from the formal condensation of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.',
        'parents': ['CHEBI:35620']
    }
}