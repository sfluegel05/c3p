"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: CHEBI:39485 ceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    Ceramides are N-acyl-sphingoid bases with a fatty acid chain (14-26 carbons) 
    attached via an amide bond to a sphingoid base containing hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find amide groups [C](=O)-[NH]
    amide_pattern = Chem.MolFromSmarts('[CX3](=O)[NX3H]')
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"

    for amide_match in amide_matches:
        # Atoms in amide_match: [C=O, N]
        carbonyl_atom = mol.GetAtomWithIdx(amide_match[0])
        nitrogen_atom = mol.GetAtomWithIdx(amide_match[1])

        # Get fatty acid chain (R-CO-)
        # Find R group attached to carbonyl (excluding the N)
        r_group = None
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetIdx() != nitrogen_atom.GetIdx() and neighbor.GetAtomicNum() == 6:
                r_group = neighbor
                break
        if not r_group:
            continue

        # Calculate fatty acid chain length (including carbonyl)
        chain_length = 1  # Start with carbonyl carbon
        visited = set([carbonyl_atom.GetIdx()])
        stack = [(r_group, 1)]  # (atom, current_length)
        max_length = 0
        while stack:
            atom, length = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6:
                max_length = max(max_length, length)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetBondTypeTo(atom) in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
                        stack.append((neighbor, length + 1))
        total_length = max_length + 1  # +1 for carbonyl
        if not (14 <= total_length <= 26):
            continue

        # Check sphingoid base (attached to N)
        # Find the carbon attached to N (excluding carbonyl)
        sphingoid_start = None
        for neighbor in nitrogen_atom.GetNeighbors():
            if neighbor.GetIdx() != carbonyl_atom.GetIdx():
                sphingoid_start = neighbor
                break
        if not sphingoid_start:
            continue

        # Check for hydroxyl groups in sphingoid base
        hydroxyl_count = 0
        visited_sph = set()
        stack_sph = [sphingoid_start]
        while stack_sph:
            atom = stack_sph.pop()
            if atom.GetIdx() in visited_sph:
                continue
            visited_sph.add(atom.GetIdx())
            # Check for hydroxyl (-OH)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                    hydroxyl_count += 1
            # Add neighboring carbons
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited_sph:
                    stack_sph.append(neighbor)

        if hydroxyl_count >= 1:
            # Check sphingoid chain length (at least 12 carbons)
            sph_chain_length = 0
            visited_chain = set()
            stack_chain = [sphingoid_start]
            while stack_chain:
                atom = stack_chain.pop()
                if atom.GetIdx() in visited_chain:
                    continue
                visited_chain.add(atom.GetIdx())
                if atom.GetAtomicNum() == 6:
                    sph_chain_length += 1
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 6 and neighbor.GetBondTypeTo(atom) in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
                            stack_chain.append(neighbor)
            if sph_chain_length >= 12:
                return True, f"Contains N-acyl (C{total_length}) sphingoid base with {hydroxyl_count} hydroxyls"

    return False, "Does not meet ceramide criteria"