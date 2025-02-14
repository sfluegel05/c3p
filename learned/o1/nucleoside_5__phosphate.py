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

    # General purine base pattern (allowing for substitutions)
    purine_smarts = 'c1ncnc2nccc1-2'
    purine_pattern = Chem.MolFromSmarts(purine_smarts)

    # General pyrimidine base pattern (allowing for substitutions)
    pyrimidine_smarts = 'c1cncnc1'
    pyrimidine_pattern = Chem.MolFromSmarts(pyrimidine_smarts)

    # Sugar (ribose or deoxyribose) pattern with possible modifications
    sugar_smarts = '[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O'
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)

    # Phosphate group pattern (allowing for mono-, di-, tri-, or tetra-phosphorylated)
    phosphate_smarts = 'OP(=O)(O)[O,P]'
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)

    # Validate patterns
    if None in [purine_pattern, pyrimidine_pattern, sugar_pattern, phosphate_pattern]:
        return False, "Error in SMARTS pattern definitions"

    # Check for purine or pyrimidine base
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)

    if not (has_purine or has_pyrimidine):
        return False, "No purine or pyrimidine base found"

    # Check for sugar ring
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No ribose or deoxyribose sugar found"

    # Check for phosphate group attached to sugar
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Verify connection between base and sugar via N-glycosidic bond
    base_atoms = set()
    if has_purine:
        base_atoms = set(mol.GetSubstructMatch(purine_pattern))
    elif has_pyrimidine:
        base_atoms = set(mol.GetSubstructMatch(pyrimidine_pattern))

    connected = False
    for sugar_match in sugar_matches:
        sugar_atoms = set(sugar_match)
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if ((begin_idx in base_atoms and end_idx in sugar_atoms) or
                (end_idx in base_atoms and begin_idx in sugar_atoms)):
                # Check if bond is between nitrogen (base) and carbon (sugar)
                atom1 = mol.GetAtomWithIdx(begin_idx)
                atom2 = mol.GetAtomWithIdx(end_idx)
                if (atom1.GetAtomicNum() == 7 and atom2.GetAtomicNum() == 6) or \
                   (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 7):
                    connected = True
                    break
        if connected:
            break

    if not connected:
        return False, "Base is not connected to sugar via N-glycosidic bond"

    # Verify phosphate groups are attached to 5' carbon of sugar
    phosphate_connected = False
    for phosphate_match in phosphate_matches:
        phosphate_atoms = set(phosphate_match)
        for sugar_match in sugar_matches:
            sugar_atoms = set(sugar_match)
            # Identify 5' carbon in sugar (has a hydroxymethyl group)
            for atom_idx in sugar_atoms:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:
                    # Look for phosphate connected to this carbon
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() in phosphate_atoms:
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