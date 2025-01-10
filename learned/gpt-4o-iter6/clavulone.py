"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    A clavulone is a class of esterified prostanoids obtained from marine corals.
    
    It's characterized by complex cyclic structures, multiple ester groups, conjugated diene
    arrangements, and may contain halogens.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General ester pattern more flexible
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Insufficient ester groups found"
    
    # Define more flexible cyclic arrangements potentially combined with esters
    cyclic_patterns = [
        Chem.MolFromSmarts("C1=CC=CCC1"),  # Simple aromatic cycles
        Chem.MolFromSmarts("C1CCC(C1)=O"), # Cyclohexanone as representative
    ]
    has_cycle = any(mol.HasSubstructMatch(cp) for cp in cyclic_patterns)
    if not has_cycle:
        return False, "No typical cyclic structure found"

    # Check for conjugated diene flexibility within the cycles
    diene_patterns = [
        Chem.MolFromSmarts("C=CC=C"),      # Conjugated diene
        Chem.MolFromSmarts("C=CCC=C"),     # Extended conjugation
    ]
    has_diene = any(mol.HasSubstructMatch(dp) for dp in diene_patterns)
    if not has_diene:
        return False, "No suitable conjugated diene arrangement found"

    # Recognize halogens presence indicative of specific clavulone subtypes
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I]")
    has_halogen = mol.HasSubstructMatch(halogen_pattern)
    
    if has_halogen:
        return True, "Matches clavulone structure with halogens"
    else:
        return True, "Matches clavulone structure possibly without halogens"

# Improved with nuances for each specification for clavulone structures.