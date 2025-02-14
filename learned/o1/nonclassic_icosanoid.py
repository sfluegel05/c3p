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

    # Count number of carbon atoms (flexible range around 20)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not 18 <= c_count <= 24:
        return False, f"Molecule has {c_count} carbons, expected between 18 and 24"

    # Count number of oxygen atoms (must have oxygenated functional groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, f"Molecule has {o_count} oxygen atoms, expected at least 2"

    # Check for long hydrocarbon chain (indicative of fatty acid backbone)
    chain_lengths = []
    paths = Chem.rdmolops.FindAllPathsOfLengthN(mol, 15, useBonds=False)
    if not paths:
        return False, "No long hydrocarbon chain of at least 15 atoms found"

    # Count number of double bonds (polyunsaturation)
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if num_double_bonds < 3:
        return False, f"Molecule has {num_double_bonds} double bonds, expected at least 3"

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

    # Ketone groups (C=O)
    ketone_pattern = Chem.MolFromSmarts('C(=O)[C!O]')
    num_ketones = len(mol.GetSubstructMatches(ketone_pattern))
    functional_groups += num_ketones

    # Peroxy groups (-OOH)
    peroxy_pattern = Chem.MolFromSmarts('OO')
    num_peroxides = len(mol.GetSubstructMatches(peroxy_pattern))
    functional_groups += num_peroxides

    if functional_groups < 2:
        return False, f"Molecule has {functional_groups} oxygenated functional groups, expected at least 2"

    # Exclude prostanoids (prostaglandins) by identifying specific cyclopentane ring
    prostaglandin_pattern = Chem.MolFromSmarts('C1C=CC(C1)CCCC(=O)O')
    if mol.HasSubstructMatch(prostaglandin_pattern):
        return False, "Molecule contains prostaglandin-like cyclopentane ring"

    # Exclude leukotrienes by identifying their conjugated triene structure
    leukotriene_pattern = Chem.MolFromSmarts('C=C/C=C/C=C')
    triene_matches = mol.GetSubstructMatches(leukotriene_pattern)
    # Check if triene is part of a specific leukotriene backbone
    for match in triene_matches:
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        if all(atom.GetDegree() == 2 for atom in atoms):  # Linear conjugated triene
            return False, "Molecule contains leukotriene-like conjugated triene"

    # Accept molecules with terminal carboxylic acids or their derivatives
    terminal_groups = [
        Chem.MolFromSmarts('C(=O)[O;H1]'),   # Carboxylic acid
        Chem.MolFromSmarts('C(=O)O[CH3]'),   # Methyl ester
        Chem.MolFromSmarts('C(=O)O[C]'),     # Ester
        Chem.MolFromSmarts('C(=O)N'),        # Amide
    ]
    has_terminal_group = any(mol.HasSubstructMatch(patt) for patt in terminal_groups)
    if not has_terminal_group:
        return False, "No terminal carboxylic acid or derivative group found"

    # If molecule passes all checks, classify as nonclassic icosanoid
    return True, "Molecule meets criteria for nonclassic icosanoid"