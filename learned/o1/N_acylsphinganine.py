"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: CHEBI:52634 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Define sphinganine backbone pattern:
    # Long-chain amino diol with hydroxyls at positions 1 and 3, amino group at position 2

    # SMARTS pattern for sphinganine backbone
    sphinganine_smarts = '[#6]-[CH](O)-[CH](N)-[CH](O)-[CH2]-[*]'
    sphinganine_pattern = Chem.MolFromSmarts(sphinganine_smarts)
    if sphinganine_pattern is None:
        return False, "Invalid sphinganine SMARTS pattern"

    # Check for sphinganine backbone
    sphinganine_matches = mol.GetSubstructMatches(sphinganine_pattern)
    if not sphinganine_matches:
        return False, "No sphinganine backbone found"

    # Define N-acyl amide group attached to the amino group
    n_acyl_smarts = '[N;H0;D2]-C(=O)-[C;!R]'
    n_acyl_pattern = Chem.MolFromSmarts(n_acyl_smarts)
    if n_acyl_pattern is None:
        return False, "Invalid N-acyl SMARTS pattern"

    # Check for N-acyl group connected to the amino nitrogen in sphinganine
    acyl_found = False
    for match in sphinganine_matches:
        amino_nitrogen_idx = match[2]
        amino_nitrogen = mol.GetAtomWithIdx(amino_nitrogen_idx)
        # Check if amino nitrogen is part of an amide
        for neighbor in amino_nitrogen.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(amino_nitrogen_idx, neighbor.GetIdx())
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and neighbor.GetAtomicNum() == 6:
                # Check if this carbon is a carbonyl carbon
                for nbr in neighbor.GetNeighbors():
                    if nbr.GetIdx() == amino_nitrogen_idx:
                        continue
                    if nbr.GetAtomicNum() == 8:
                        bond2 = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nbr.GetIdx())
                        if bond2.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            # Found amide bond
                            # Now check the length of the acyl chain
                            acyl_chain_length = 0
                            visited = set()
                            queue = [nbr for nbr in neighbor.GetNeighbors() if nbr.GetIdx() != amino_nitrogen_idx and nbr.GetAtomicNum() == 6]
                            while queue:
                                atom = queue.pop(0)
                                atom_idx = atom.GetIdx()
                                if atom_idx in visited:
                                    continue
                                visited.add(atom_idx)
                                acyl_chain_length += 1
                                for nbr in atom.GetNeighbors():
                                    if nbr.GetIdx() not in visited and nbr.GetAtomicNum() == 6 and not nbr.IsInRing():
                                        queue.append(nbr)
                            if acyl_chain_length >= 12:
                                acyl_found = True
                                break
            if acyl_found:
                break
        if acyl_found:
            break

    if not acyl_found:
        return False, "No N-acyl group with sufficient length found attached to amino nitrogen"

    return True, "Contains sphinganine backbone with N-acyl group forming an amide bond"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:52634',
        'name': 'N-acylsphinganine',
        'definition': "A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.",
        'parents': ['CHEBI:33598', 'CHEBI:33709']
    },
    'config': {}
}