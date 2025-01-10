"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is any fatty acid ester where the carboxylic acid component is lauric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the lauric acid moiety pattern (C12 chain with a carboxyl group)
    lauric_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC(=O)[OX2]")
    if not mol.HasSubstructMatch(lauric_acid_pattern):
        return False, "No lauric acid moiety found"

    # Check for ester bond (-COO-)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester bond found"

    # Verify that the ester bond is connected to the lauric acid moiety
    lauric_acid_match = mol.GetSubstructMatch(lauric_acid_pattern)
    lauric_acid_atoms = set(lauric_acid_match)
    ester_atoms = set()
    for match in ester_matches:
        ester_atoms.update(match)
    
    if not lauric_acid_atoms.intersection(ester_atoms):
        return False, "Ester bond not connected to lauric acid moiety"

    return True, "Contains lauric acid moiety with an ester bond"