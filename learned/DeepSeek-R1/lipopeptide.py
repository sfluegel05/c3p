"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: CHEBI: ??? lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide consists of a peptide (multiple amide bonds) with an attached lipid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for peptide component: at least two amide bonds
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 2:
        return False, f"Found {len(amide_matches)} amide groups, need at least 2 for a peptide"

    # Precompute ring information to exclude cyclic chains
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Function to calculate maximum linear aliphatic carbon chain length from a given atom
    def max_carbon_chain(atom, visited=None):
        if visited is None:
            visited = set()
        # Skip if atom is in a ring or already visited
        if atom.GetAtomicNum() != 6 or atom.GetIdx() in visited:
            return 0
        for ring in rings:
            if atom.GetIdx() in ring:
                return 0
        visited.add(atom.GetIdx())
        max_length = 0
        for neighbor in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond.GetBondType() == Chem.BondType.SINGLE and not bond.IsInRing():
                length = 1 + max_carbon_chain(neighbor, visited)
                if length > max_length:
                    max_length = length
        visited.remove(atom.GetIdx())
        return max_length

    # Check amide groups for lipid chains (minimum 6 carbons)
    lipid_found = False
    for amide_match in amide_matches:
        n_idx = amide_match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        # Check substituents on the amide nitrogen (excluding the carbonyl carbon)
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                chain_length = max_carbon_chain(neighbor)
                if chain_length >= 5:  # 5 carbons from neighbor = total 6
                    lipid_found = True
                    break
        if lipid_found:
            break

    # Check ester groups for lipid chains
    if not lipid_found:
        ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        for ester_match in ester_matches:
            o_idx = ester_match[0]
            o_atom = mol.GetAtomWithIdx(o_idx)
            # Check substituents on the ester oxygen (excluding the carbonyl carbon)
            for neighbor in o_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    chain_length = max_carbon_chain(neighbor)
                    if chain_length >= 5:
                        lipid_found = True
                        break
            if lipid_found:
                break

    # Check ether groups for lipid chains
    if not lipid_found:
        ether_pattern = Chem.MolFromSmarts("[OX2]")
        ether_matches = mol.GetSubstructMatches(ether_pattern)
        for ether_match in ether_matches:
            o_idx = ether_match[0]
            o_atom = mol.GetAtomWithIdx(o_idx)
            # Check both sides of the ether oxygen
            for neighbor in o_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    chain_length = max_carbon_chain(neighbor)
                    if chain_length >= 5:
                        lipid_found = True
                        break
            if lipid_found:
                break

    if not lipid_found:
        return False, "No lipid chain (>=6 carbons) attached via amide, ester, or ether"

    return True, "Contains peptide with lipid chain attached via amide/ester/ether"