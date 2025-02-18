"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are oxidation products of C20 essential fatty acids (EFAs) like arachidonic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check total carbons >=18 (allowing for some modification from C20 EFAs)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:
        return False, f"Only {c_count} carbons, need at least 18"

    # Check for carboxylic acid (-COOH) or ester (-COO-) groups
    carboxylic_acid = mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3](=O)[OX2H1]'))
    ester = mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3](=O)[OX2][#6]'))
    if not (carboxylic_acid or ester):
        return False, "No carboxylic acid or ester group"

    # Check for at least two double bonds (common in EFA derivatives)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 2:
        return False, f"Only {double_bonds} double bonds, need at least 2"

    # Check for oxygen-containing functional groups (hydroxyl, epoxide, etc.)
    hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H]'))
    epoxide = mol.HasSubstructMatch(Chem.MolFromSmarts('[O]1CCO1'))
    peroxide = mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2][OX2]'))
    ketone = mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3]=O'))
    ether = mol.HasSubstructMatch(Chem.MolFromSmarts('[OD2]([#6])[#6]'))
    
    if not any([hydroxyl, epoxide, peroxide, ketone, ether]):
        return False, "No oxygen-containing groups (hydroxyl/epoxide/peroxide/ketone/ether)"

    # Check total oxygen atoms >=2 (oxidation products)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, f"Only {o_count} oxygen atoms, need at least 2"

    return True, "Derived from C20 EFAs with oxidation markers (long chain, oxygen groups, multiple double bonds)"