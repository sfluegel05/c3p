"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI fatty amide (monocarboxylic acid amide derived from a fatty acid)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid (aliphatic chain >=4 carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find all amide groups (including primary/secondary/tertiary)
    amide_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[NX3H2,H1,H0]')
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"

    # Check for carboxylic acid groups (exclude if present)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Contains carboxylic acid group"

    # Helper function to calculate aliphatic chain length via BFS
    def get_aliphatic_chain_length(start_atom):
        visited = set()
        max_length = 0
        queue = [(start_atom, 1)]  # (atom, current_length)
        while queue:
            atom, length = queue.pop(0)
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            # Skip non-carbon, aromatic, or ring atoms
            if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.IsInRing():
                continue
            if length > max_length:
                max_length = length
            # Traverse single bonds to aliphatic carbons
            for neighbor in atom.GetNeighbors():
                bond = atom.GetBondToAtom(neighbor)
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    queue.append((neighbor, length + 1))
        return max_length

    # Check each amide group for qualifying chain
    for match in amide_matches:
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Find R-group carbon attached to carbonyl
        r_carbon = None
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in {match[1], match[2]}:
                r_carbon = neighbor
                break
        if not r_carbon:
            continue  # No R-group found for this amide

        chain_length = get_aliphatic_chain_length(r_carbon)
        if chain_length >= 4:
            return True, f"Amide with {chain_length}-carbon aliphatic chain"

    return False, "No amide with sufficient aliphatic chain length"