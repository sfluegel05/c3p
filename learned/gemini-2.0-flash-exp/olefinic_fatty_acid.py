"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Any fatty acid containing at least one C=C double bond.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is a fatty acid with at least one C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    if not _check_carboxylic_acid(mol):
        return False, "No carboxylic acid group found"
    
    if not _check_double_bond(mol):
         return False, "No C=C double bond found"

    if not _check_fatty_acid_properties(mol):
        return False, "Molecule does not meet fatty acid property requirements"

    return True, "Contains a carboxylic acid group and at least one C=C double bond and has 8 or more carbons, enough rotatable bonds and a reasonable weight."

def _check_carboxylic_acid(mol):
    """
    Checks if the molecule contains a carboxylic acid group, accounting for implicit hydrogens
    """
    acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if mol.HasSubstructMatch(acid_pattern):
        return True
    
    #check for O=C-O without H
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    for neighbor2 in atom.GetNeighbors():
                        if neighbor2.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom.GetIdx(),neighbor2.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                            return True
    return False


def _check_double_bond(mol):
    """
    Checks if the molecule contains at least one C=C double bond.
    """
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    return mol.HasSubstructMatch(double_bond_pattern)


def _check_fatty_acid_properties(mol):
    """
    Checks molecule properties to ensure it looks like a fatty acid, such as carbon count, rotatable bonds, weight, and oxygen count.
    """
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 8:
         return False
    
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rotatable_bonds < 4:
        return False
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False
    
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 4:
        return False
    
    return True