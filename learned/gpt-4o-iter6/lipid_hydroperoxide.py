"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is any lipid carrying one or more hydroperoxy (-OOH) substituents.

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

    # Look for hydroperoxy groups (-OOH)
    hydroperoxy_pattern = Chem.MolFromSmarts("[O][O][H]")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hydroperoxy_matches:
        return False, "No hydroperoxy group found"
    
    # Improve lipid detection by refining carbon chain pattern and checking for common lipid features
    # Check for long carbon chains including possible double bonds and branching
    carbon_chain_pattern = Chem.MolFromSmarts("[C](=[C])[C]")  # Assumes presence of bonds typical of lipids
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if carbon_chain_matches:
        # Checking for long carbon chain with bonds resembling lipid characteristics
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
        
        if carbon_count >= 10 and double_bond_count >= 2:
            return True, "Contains hydroperoxy group(s) within a lipid structure"
        
    # Checking for the presence of end groups like carboxylic acid that are common in lipids
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return True, "Contains hydroperoxy group(s) within a lipid structure"
    
    return False, "Does not fit lipid hydroperoxide profile"