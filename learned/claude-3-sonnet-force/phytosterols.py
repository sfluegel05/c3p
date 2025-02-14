"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:36177 phytosterol

Phytosterols are sterols similar to cholesterol which occur in plants and vary only
in carbon side chains and/or presence or absence of a double bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phytosterol(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("[C@H]1[C@H]2[C@@H]([C@]([C@@H]([C@H]([C@@H]2[C@H]([C@@H]1O)C)C)C)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for number of carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if not (27 <= c_count <= 32 and o_count == 1):
        return False, f"Incorrect number of carbons ({c_count}) or oxygens ({o_count}) for phytosterol"
    
    # Check for double bonds and aliphatic side chains
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    if double_bonds == 0 or rotatable_bonds < 3:
        return False, "Missing double bonds or aliphatic side chains"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (350 <= mol_wt <= 550):
        return False, f"Molecular weight {mol_wt} outside typical range for phytosterols"
    
    return True, "Contains steroid backbone with aliphatic side chains and double bonds"