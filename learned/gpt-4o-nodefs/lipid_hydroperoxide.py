"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for hydroperoxide group (-OOH)
    hydroperoxide_pattern = Chem.MolFromSmarts("[OH,OX2]-[OX1H0,R0]")
    if not mol.HasSubstructMatch(hydroperoxide_pattern):
        return False, "No hydroperoxide group found"
    
    # Check for a carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for long carbon chains (at least 12 carbons)
    carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_chain_length < 12:
        return False, f"Carbon chain is too short, found {carbon_chain_length} carbons"
    
    # Check for multiple double bonds (polyunsaturation)
    double_bonds = len([bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE])
    if double_bonds < 2:
        return False, f"Insufficient double bonds, found {double_bonds}"
    
    return True, "Molecule matches structure of a lipid hydroperoxide"