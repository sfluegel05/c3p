"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool or None: True if molecule is a non-proteinogenic amino acid, False if not, None if unable to determine
        str or None: Reason for classification, or None if unable to determine
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for the basic amino acid backbone: an amino (-NH2) group and a carboxylic acid (-COOH) group
    amino_group_pattern = Chem.MolFromSmarts("[NX3][CX4]")
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OX1H1]")
    
    if not mol.HasSubstructMatch(amino_group_pattern) or not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Lacks standard amino acid backbone elements"

    # Define patterns for unusual side chains or modifications that indicate non-proteinogenic status
    uncommon_side_chain_patterns = [
        Chem.MolFromSmarts("[CX4]([#6])([#6])([#9]),[OX2][CX3](=[OX1])"), # Fluorinated groups, esters
        Chem.MolFromSmarts("[CX4]#N"),  # Cyanide group
        Chem.MolFromSmarts("[Se]"),  # Selenium-containing
        Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX3]=[NX2]"), # Imines or other unusual nitrogen groups
        # Add other patterns indicative of non-standard amino acids
    ]

    for pattern in uncommon_side_chain_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains unusual modification or uncommon side chain"

    # Additional checks for modified backbone or unique configurations could go here

    # If none of the patterns match, fall back to indicating uncertainty
    return None, "Could not determine if it is a non-proteinogenic amino acid"