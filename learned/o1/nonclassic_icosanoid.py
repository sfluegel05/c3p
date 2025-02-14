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

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Allow molecules with 20 to 25 carbons to account for side chains
    if not 20 <= c_count <= 25:
        return False, f"Molecule has {c_count} carbons, expected between 20 and 25"

    # Check for terminal carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No terminal carboxylic acid group found"

    # Exclude molecules with peptide bonds (amide linkages between carbons)
    peptide_bond = Chem.MolFromSmarts('C(=O)N[C]')
    if mol.HasSubstructMatch(peptide_bond):
        return False, "Molecule contains peptide bonds (possible peptide/protein)"

    # Exclude prostaglandins by identifying cyclopentane ring with side chains
    prostaglandin_pattern = Chem.MolFromSmarts('C1CCCC1')
    if mol.HasSubstructMatch(prostaglandin_pattern):
        # Further check for side chains to confirm prostaglandin
        side_chain = Chem.MolFromSmarts('C1CCCC1CC')
        if mol.HasSubstructMatch(side_chain):
            return False, "Molecule contains prostaglandin-like cyclopentane ring"

    # Count total number of oxygen atoms
    o_count_total = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    # Count number of carboxylic acid groups
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid)
    o_count_carboxy = len(carboxy_matches) * 2  # Each carboxylic acid has 2 oxygens

    # Calculate number of oxygens beyond carboxylic acid(s)
    o_count_additional = o_count_total - o_count_carboxy
    if o_count_additional < 1:
        return False, "No additional oxygenated functional groups found"

    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))

    # Check for epoxides (three-membered cyclic ethers)
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    num_epoxides = len(mol.GetSubstructMatches(epoxide_pattern))

    # Exclude molecules with ketone groups (C=O not part of carboxylic acid)
    ketone_pattern = Chem.MolFromSmarts('[CX3](=O)[#6]')
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    # Exclude ketones that are not part of carboxylic acids
    num_ketones = 0
    for match in ketone_matches:
        ketone_atom = mol.GetAtomWithIdx(match[0])
        is_in_carboxy = any(ketone_atom.IsInRing() for match in carboxy_matches)
        if not is_in_carboxy:
            num_ketones += 1
    if num_ketones > 0:
        return False, "Molecule contains ketone groups (possible prostanoid)"

    # Sum up oxygenated functional groups
    num_oxygenated_groups = num_hydroxyls + num_epoxides
    if num_oxygenated_groups < 1:
        return False, "No significant oxygenated functional groups found"

    # Exclude molecules with multiple rings unrelated to epoxides
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > num_epoxides:
        return False, "Molecule contains rings unrelated to epoxides"

    # Check for long hydrocarbon chain
    # Identify the longest carbon chain in the molecule
    chains = Chem.GetSSSR(mol)
    max_chain_length = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                # Simple heuristic to estimate chain length
                chain_length = len(Chem.rdmolops.GetShortestPath(mol, atom1.GetIdx(), atom2.GetIdx()))
                if chain_length > max_chain_length:
                    max_chain_length = chain_length
    if max_chain_length < 10:
        return False, "No long hydrocarbon chain found"

    # If molecule passes all checks, classify as nonclassic icosanoid
    return True, "Molecule meets criteria for nonclassic icosanoid"