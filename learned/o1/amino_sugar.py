"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:17754 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is any sugar having one or more alcoholic hydroxy groups
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for monosaccharides (both cyclic and acyclic forms)
    # Pattern for monosaccharide rings (furanose and pyranose)
    sugar_ring_smarts = "[C,O]1[C,O][C,O][C,O][C,O][C,O]1"  # Pyranose (6-membered ring)
    sugar_ring_smarts_5 = "[C,O]1[C,O][C,O][C,O][C,O]1"     # Furanose (5-membered ring)
    sugar_ring_pattern = Chem.MolFromSmarts(sugar_ring_smarts)
    sugar_ring_pattern_5 = Chem.MolFromSmarts(sugar_ring_smarts_5)

    # Pattern for anomeric carbon with amino group replacing OH
    amino_substitution_smarts = "[C;R][N;D2]"
    amino_substitution_pattern = Chem.MolFromSmarts(amino_substitution_smarts)

    # Flag to track if amino sugar is found
    is_amino_sugar_found = False

    # Check for cyclic amino sugars
    # Find all sugar rings
    ring_matches = mol.GetSubstructMatches(sugar_ring_pattern) + mol.GetSubstructMatches(sugar_ring_pattern_5)
    if ring_matches:
        for ring in ring_matches:
            ring_atoms = list(ring)
            ring_atom_indices = set(ring_atoms)
            # Check each carbon in the ring for amino substitution
            for idx in ring_atoms:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # Carbon atom
                    # Check for original hydroxy group (carbon attached to oxygen)
                    has_oxygen = False
                    has_amino = False
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8:
                            has_oxygen = True
                        elif neighbor.GetAtomicNum() == 7:
                            has_amino = True
                    # If amino group is present and replaces hydroxy group
                    if has_amino:
                        is_amino_sugar_found = True
                        return True, "Contains sugar ring with amino group(s) replacing hydroxy group(s)"

    # Check for acyclic amino sugars
    # Acyclic monosaccharide pattern: chain of 4-6 carbons with hydroxy groups
    acyclic_sugar_smarts = "[CX4H1;!$(C=*)][CX4H1][CX4H1][CX4H1][CX4H1][CX4H1]"
    acyclic_sugar_pattern = Chem.MolFromSmarts(acyclic_sugar_smarts)
    acyclic_matches = mol.GetSubstructMatches(acyclic_sugar_pattern)
    if not acyclic_matches:
        # Try shorter chains
        acyclic_sugar_smarts = "[CX4H1;!$(C=*)][CX4H1][CX4H1][CX4H1][CX4H1]"
        acyclic_sugar_pattern = Chem.MolFromSmarts(acyclic_sugar_smarts)
        acyclic_matches = mol.GetSubstructMatches(acyclic_sugar_pattern)

    if acyclic_matches:
        for match in acyclic_matches:
            match_atoms = list(match)
            amino_found = False
            for idx in match_atoms:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:
                    has_hydroxy = False
                    has_amino = False
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8:
                            # Check if it's a hydroxyl group
                            if neighbor.GetTotalDegree() == 1:
                                has_hydroxy = True
                        elif neighbor.GetAtomicNum() == 7:
                            has_amino = True
                    # If amino group is present and replaces hydroxy group
                    if has_amino:
                        amino_found = True
            if amino_found:
                is_amino_sugar_found = True
                return True, "Contains acyclic amino sugar with amino group(s) replacing hydroxy group(s)"

    # Additional check: amino sugar with open-chain aldehyde or ketone
    open_chain_sugar_smarts = "[C;X4][C;X4][C;X4][C;X4][C;X4]"
    open_chain_sugar_pattern = Chem.MolFromSmarts(open_chain_sugar_smarts)
    open_chain_matches = mol.GetSubstructMatches(open_chain_sugar_pattern)
    if open_chain_matches:
        for match in open_chain_matches:
            amino_found = False
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:
                    has_hydroxy = False
                    has_amino = False
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8:
                            # Check for hydroxyl group
                            if neighbor.GetTotalDegree() == 1:
                                has_hydroxy = True
                        elif neighbor.GetAtomicNum() == 7:
                            has_amino = True
                    if has_amino:
                        amino_found = True
            if amino_found:
                is_amino_sugar_found = True
                return True, "Contains open-chain amino sugar with amino group(s) replacing hydroxy group(s)"

    if is_amino_sugar_found:
        return True, "Contains amino sugar structure"

    return False, "Does not contain amino sugar ring or acyclic amino sugar"