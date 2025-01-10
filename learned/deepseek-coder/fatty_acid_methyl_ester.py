"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: CHEBI:59826 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is a carboxylic ester obtained by the formal condensation of a fatty acid with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the ester group pattern: [CX3](=[OX1])[OX2][CH3]
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) < 1:
        return False, "No methyl ester group found"

    # Check for a fatty acid chain (long carbon chain attached to the ester)
    # We look for a carbon chain with at least 6 carbons
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    
    if len(fatty_acid_matches) < 1:
        return False, "No fatty acid chain found"

    # Count the number of carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if c_count < 8:
        return False, "Too few carbons to be a fatty acid methyl ester"

    # Check molecular weight - fatty acid methyl esters typically have MW > 150 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for a fatty acid methyl ester"

    # Additional check to ensure the ester is part of a carboxylic acid ester
    carboxylic_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3]")
    carboxylic_ester_matches = mol.GetSubstructMatches(carboxylic_ester_pattern)
    
    if len(carboxylic_ester_matches) < 1:
        return False, "Ester group is not part of a carboxylic acid ester"

    # Ensure the ester is connected to a long carbon chain
    ester_connected_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX3](=[OX1])[OX2][CH3]")
    ester_connected_matches = mol.GetSubstructMatches(ester_connected_pattern)
    
    if len(ester_connected_matches) < 1:
        return False, "Ester group is not connected to a long carbon chain"

    return True, "Contains a methyl ester group attached to a fatty acid chain"