"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine
    or purine base in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define generalized SMARTS patterns

    # Generalized purine or substituted purine base pattern
    purine_smarts = 'c1ncnc2ncccc12'
    purine_pattern = Chem.MolFromSmarts(purine_smarts)

    # Generalized pyrimidine or substituted pyrimidine base pattern
    pyrimidine_smarts = 'c1cncnc1'
    pyrimidine_pattern = Chem.MolFromSmarts(pyrimidine_smarts)

    # Generalized sugar ring pattern (furanose ring)
    sugar_smarts = '[C,R]1([O,R])[C,R][C,R][C,R][O,R]1'
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)

    # Phosphate group pattern (mono-, di-, tri-, or tetra-phosphate)
    phosphate_smarts = 'OP(=O)([O,P])([O,P])[O,P]'
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)

    # Validate patterns
    if None in [purine_pattern, pyrimidine_pattern, sugar_pattern, phosphate_pattern]:
        return False, "Error in SMARTS pattern definitions"

    # Check for purine or pyrimidine base (including substitutions)
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)

    if not (has_purine or has_pyrimidine):
        return False, "No purine or pyrimidine base found"

    # Check for sugar ring (ribose or deoxyribose including modifications)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar ring (furanose) found"

    # Check for phosphate group attached to sugar
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Verify connection between base and sugar via glycosidic bond
    base_atoms = set()
    if has_purine:
        base_atoms = set(mol.GetSubstructMatch(purine_pattern))
    elif has_pyrimidine:
        base_atoms = set(mol.GetSubstructMatch(pyrimidine_pattern))

    connected = False
    for sugar_match in sugar_matches:
        sugar_atoms = set(sugar_match)
        # Check for bond between base and sugar
        for atom_idx in base_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in sugar_atoms:
                    connected = True
                    break
            if connected:
                break
        if connected:
            break

    if not connected:
        return False, "Base is not connected to sugar"

    # Verify phosphate group is attached to 5' carbon of sugar
    phosphate_connected = False
    for sugar_match in sugar_matches:
        sugar_atoms = set(sugar_match)
        # Identify 5' carbon in sugar ring
        five_prime_carbons = []
        for atom_idx in sugar_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                # Look for hydroxymethyl group attached to sugar ring
                num_ring_bonds = sum(1 for bond in atom.GetBonds() if bond.IsInRing())
                if num_ring_bonds == 1:
                    five_prime_carbons.append(atom)
        # Check if phosphate is connected to any of the 5' carbons
        for five_c_atom in five_prime_carbons:
            for neighbor in five_c_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    for neighbor2 in neighbor.GetNeighbors():
                        if neighbor2.GetIdx() in phosphate_matches[0]:
                            phosphate_connected = True
                            break
                    if phosphate_connected:
                        break
            if phosphate_connected:
                break
        if phosphate_connected:
            break

    if not phosphate_connected:
        return False, "Phosphate group is not attached to 5' carbon of sugar"

    # At this point, criteria are satisfied
    return True, "Molecule is a nucleoside 5'-phosphate"