"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    A nonclassic icosanoid is a biologically active signaling molecule made by oxygenation
    of C20 fatty acids other than the classic icosanoids (the leukotrienes and the prostanoids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for atoms other than C, H, and O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, f"Molecule contains heteroatom {atom.GetSymbol()}, expected only C, H, and O"

    # Exclude molecules with aromatic rings
    if mol.GetRingInfo().IsAtomInRingOfSize(0, 6):
        aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
        if aromatic_atoms:
            return False, "Molecule contains aromatic rings"

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not 16 <= c_count <= 26:
        return False, f"Molecule has {c_count} carbons, expected between 16 and 26"

    # Check for long aliphatic chain (indicative of fatty acid backbone)
    chains = Chem.GetSymmSSSR(mol)
    longest_chain = 0
    for path in Chem.rdmolops.FindAllPathsOfLengthN(mol, 15, useBonds=False):
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path):
            chain_length = len(path)
            if chain_length > longest_chain:
                longest_chain = chain_length
    if longest_chain < 15:
        return False, "No long hydrocarbon chain of at least 15 carbons found"

    # Count number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, f"Molecule has {o_count} oxygen atoms, expected at least 2"

    # Check for oxygenated functional groups
    functional_groups = 0

    # Hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    functional_groups += num_hydroxyls

    # Epoxide rings (three-membered cyclic ethers)
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    num_epoxides = len(mol.GetSubstructMatches(epoxide_pattern))
    functional_groups += num_epoxides

    if functional_groups < 1:
        return False, "Molecule lacks characteristic oxygenated groups (e.g., hydroxyls, epoxides)"

    # Exclude prostanoids by identifying cyclopentane ring connected to aliphatic chains
    prostanoid_pattern = Chem.MolFromSmarts('C1CCCC1')
    matches = mol.GetSubstructMatches(prostanoid_pattern)
    for match in matches:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        if all(atom.GetAtomicNum() == 6 for atom in ring_atoms):
            return False, "Molecule contains cyclopentane ring characteristic of prostanoids"

    # Exclude molecules with more than one ring (excluding epoxides)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    num_epoxides = len(mol.GetSubstructMatches(epoxide_pattern))
    if num_rings - num_epoxides > 0:
        return False, "Molecule contains rings other than epoxides"

    # Check for terminal carboxylic acid or derivative
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C]')
    terminal_group_found = False
    for pattern in [carboxylic_pattern, ester_pattern]:
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            atom_idx = match[0]
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check if it's at the end of a chain
            if atom.GetDegree() == 2:
                terminal_group_found = True
                break
        if terminal_group_found:
            break
    if not terminal_group_found:
        return False, "No terminal carboxylic acid or ester group found"

    return True, "Molecule meets criteria for nonclassic icosanoid"