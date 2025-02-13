"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: CHEBI:36975 Essential fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are polyunsaturated fatty acids with multiple double bonds,
    typically found in nature and required in the human diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid (-COOH) group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "Not a carboxylic acid (missing -COOH group)"

    # Check for long carbon chain (>10 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Carbon chain too short for fatty acid"

    # Check for multiple double bonds (>=2)
    double_bond_pattern = Chem.MolFromSmarts("=C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 2:
        return False, "Not polyunsaturated (needs at least 2 double bonds)"

    # Check for cis/trans configuration
    cis_pattern = Chem.MolFromSmarts("/C=C/")
    trans_pattern = Chem.MolFromSmarts("\\C=C\\")
    if not (mol.HasSubstructMatch(cis_pattern) or mol.HasSubstructMatch(trans_pattern)):
        return False, "Double bonds lack cis/trans configuration"

    # Check for presence in natural sources
    # This is a heuristic and not a definitive check
    natural_sources = ["fish", "plant", "algae", "nut", "seed"]
    mol_name = Chem.MolToSmiles(mol, isomericSmiles=True).split("_")[0].lower()
    if not any(source in mol_name for source in natural_sources):
        return False, "Not found in natural sources (heuristic)"

    return True, "Polyunsaturated fatty acid with multiple cis/trans double bonds, found in natural sources"