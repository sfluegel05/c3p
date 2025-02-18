"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is a compound with the general structure RC(=O)H.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldehyde functional group pattern: [CX3H1](=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    
    # Check if the molecule contains the aldehyde pattern
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Additional check to ensure it's not part of a carboxylic acid or ester
        carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
        ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CX4]")
        amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
        
        # Check if the aldehyde carbon is part of a carboxylic acid, ester, or amide
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
        for match in aldehyde_matches:
            aldehyde_carbon = match[0]
            # Check if the aldehyde carbon is part of a carboxylic acid, ester, or amide
            if (mol.HasSubstructMatch(carboxylic_acid_pattern) and 
                any(atom.GetIdx() == aldehyde_carbon for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)) or \
               (mol.HasSubstructMatch(ester_pattern) and 
                any(atom.GetIdx() == aldehyde_carbon for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)) or \
               (mol.HasSubstructMatch(amide_pattern) and 
                any(atom.GetIdx() == aldehyde_carbon for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)):
                return False, "Contains carbonyl group but is part of a carboxylic acid, ester, or amide"
        
        return True, "Contains the aldehyde functional group (RC(=O)H)"
    else:
        return False, "Does not contain the aldehyde functional group (RC(=O)H)"