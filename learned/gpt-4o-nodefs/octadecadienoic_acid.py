"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid generally has an 18-carbon backbone with exactly 2 cis/trans double bonds,
    and a terminal carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal carboxylic acid group
    # (matches COOH at the end but allows for additional branching elsewhere)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"

    # Count the number of double bonds, including stereo information
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    stereochem_checks = {
        'E_isomers': Chem.MolFromSmarts("C/C=C/C"),
        'Z_isomers': Chem.MolFromSmarts("C/C=C\\C")
    }
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    e_count = len(mol.GetSubstructMatches(stereochem_checks['E_isomers']))
    z_count = len(mol.GetSubstructMatches(stereochem_checks['Z_isomers']))

    if len(double_bond_matches) != 2:
        return False, f"Expected exactly 2 double bonds, found {len(double_bond_matches)}"
    
    # Count total number of carbon atoms in the backbone only
    backbone_matches = Chem.rdmolops.GetMolFrags(mol, asMols=True)
    longest_chain = max(backbone_matches, key=lambda frag: frag.GetNumAtoms())
    c_count = sum(1 for atom in longest_chain.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if the longest carbon backbone is 18 atoms (accounting for potential branching)
    if c_count != 18:
        return False, f"Expected a main carbon backbone of 18 atoms, found {c_count}"

    return True, "Molecule is an octadecadienoic acid with a correct carbon backbone and double bonds"