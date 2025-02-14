"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are C25 isoprenoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20 or carbon_count > 30 :
        return False, f"Number of carbons is {carbon_count}, expected around 25 for a sesterterpenoid"

    # Count methyl groups (CH3) attached to a non-hydrogen
    methyl_pattern = Chem.MolFromSmarts("[CX4H3]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    methyl_count = len(methyl_matches)
    if methyl_count < 4:
       return False, f"Too few methyl groups ({methyl_count}), expected at least 4."
    
    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("[*]=[*]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    double_bond_count = len(double_bond_matches)
    if double_bond_count < 3:
       return False, f"Too few double bonds ({double_bond_count}), expected at least 3"

    # Check number of rotatable bonds
    num_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rotatable < 5:
        return False, f"Too few rotatable bonds, expected at least 5."
    
     # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 500 :
        return False, f"Molecular weight is {mol_wt}, expected between 300 and 500 for a sesterterpenoid"

    return True, "Matches criteria for a sesterterpenoid based on number of carbons, methyl groups, unsaturation, rotatable bonds and molecular weight."