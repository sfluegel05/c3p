"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine has a glycerol backbone with two fatty acid esters and a phosphoethanolamine group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define the core phosphatidylethanolamine pattern using SMARTS
    # Pattern matches glycerol backbone with two ester groups and phosphoethanolamine group
    pe_pattern = Chem.MolFromSmarts(
        "[CH2]-[CH](-[OX2]-C(=O))-[CH2](-[OX2]-C(=O))-"  # Glycerol with two ester groups
        "O-P(=O)(-O-[CH2]-[CH2]-N([CH3])[CH3])-O"  # Phosphoethanolamine group (allowing N-methylation)
    )
    
    # Allow variations in N-methylation (0, 1, or 2 methyl groups)
    pe_pattern_variants = [
        Chem.MolFromSmarts(
            "[CH2]-[CH](-[OX2]-C(=O))-[CH2](-[OX2]-C(=O))-"
            "O-P(=O)(-O-[CH2]-[CH2]-N[CH3])-O"
        ),
        Chem.MolFromSmarts(
            "[CH2]-[CH](-[OX2]-C(=O))-[CH2](-[OX2]-C(=O))-"
            "O-P(=O)(-O-[CH2]-[CH2]-NH)-O"
        )
    ]
    
    # Check main pattern or variants
    if mol.HasSubstructMatch(pe_pattern):
        return True, "Glycerol backbone with two esters and phosphoethanolamine group"
    for variant in pe_pattern_variants:
        if mol.HasSubstructMatch(variant):
            return True, "Glycerol backbone with two esters and phosphoethanolamine group"
    
    # If no matches, check for possible stereochemistry or alternate connectivity
    # More flexible pattern allowing any glycerol connectivity
    flexible_pe_pattern = Chem.MolFromSmarts(
        "[C;H1,H2]([OX2]-C(=O))([OX2]-C(=O))-C-O-P(=O)(O-[CH2]-[CH2]-N)-O"
    )
    if mol.HasSubstructMatch(flexible_pe_pattern):
        return True, "Glycerol backbone with two esters and phosphoethanolamine group"
    
    return False, "Missing required glycerol-ester-phosphoethanolamine structure"