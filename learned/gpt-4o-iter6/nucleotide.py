"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phosphate group: O=P(O)(O)O
    phosphate_pattern = Chem.MolFromSmarts("[OX2,PX4](=O)[OX1][OX1]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Identify sugar ring with nucleoside connection
    # Assuming ribose & deoxyribose as part of nucleotides: [CH2X4][Ox2][CX4](Nuc)...
    ribose_pattern = Chem.MolFromSmarts("[C@H1]([C@H1]([C@H1]([C@H1]([OX2])[OX2])[CX4])[C@H1](O)O)O")
    
    if not mol.HasSubstructMatch(ribose_pattern):
        # If ribose is not found, check for deoxyribose (missing one OH)
        deoxyribose_pattern = Chem.MolFromSmarts("[C@H1]([C@H1]([C@H1]([C@H1]([OX2])[CX4])[OX2])[C@H1](O)O)O")
        if not mol.HasSubstructMatch(deoxyribose_pattern):
            return False, "No nucleoside sugar ring found"
    
    # Check for the connection of phosphate to 3' or 5' carbon atom
    connection_pattern_3_5 = Chem.MolFromSmarts("[C@H1]([CX4][OX2][PX4](=O)([OX1])[OX1])[OX2]")
    
    if not mol.HasSubstructMatch(connection_pattern_3_5):
        return False, "No proper phosphate connection to 3' or 5' carbon"
    
    return True, "Contains nucleoside with phosphate attached to 3' or 5' carbon"