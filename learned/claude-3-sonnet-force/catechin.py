"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:29790 catechin

Catechins are members of the class of hydroxyflavan that have a flavan-3-ol skeleton
and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if not formula.startswith("C15H14O6") and not any(formula.startswith(f"C{i}H{j}O{k}") for i in range(16, 31) for j in range(14, 31) for k in range(7, 16)):
        return False, "Molecular formula does not match catechin or its substituted derivatives"

    # Check for flavonoid backbone
    flavonoid_pattern = Chem.MolFromSmarts("c1c(O)cc(O)c2c1C(O)C(O)Cc2")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "Molecule does not contain a flavonoid backbone"

    # Check for flavan-3-ol skeleton
    flavan3ol_pattern = Chem.MolFromSmarts("[C@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1ccc(O)cc1")
    if not mol.HasSubstructMatch(flavan3ol_pattern):
        return False, "Molecule does not contain a flavan-3-ol skeleton"

    # Check for common substituents and modifications
    galloyl_group = Chem.MolFromSmarts("C(=O)c1cc(O)c(O)c(O)c1")
    ester_group = Chem.MolFromSmarts("C(=O)O")
    methoxy_group = Chem.MolFromSmarts("OC")
    hydroxyl_group = Chem.MolFromSmarts("O")

    n_galloyl = len(mol.GetSubstructMatches(galloyl_group))
    n_ester = len(mol.GetSubstructMatches(ester_group))
    n_methoxy = len(mol.GetSubstructMatches(methoxy_group))
    n_hydroxyl = len(mol.GetSubstructMatches(hydroxyl_group))

    # Additional structural constraints
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 3 or n_rings > 5:
        return False, "Number of rings does not match typical catechin structures"

    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_aromatic_rings < 2:
        return False, "Insufficient number of aromatic rings for a catechin"

    if n_hydroxyl < 3:
        return False, "Too few hydroxyl groups for a catechin"

    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 6 or oxygen_count > 12:
        return False, "Number of oxygen atoms does not match typical catechin structures"

    return True, "Molecule meets the structural requirements for a catechin or its substituted derivative"