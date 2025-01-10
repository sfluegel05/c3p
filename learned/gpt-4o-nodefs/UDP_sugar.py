"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expanded SMARTS patterns
    # Uridine includes both pyrimidine and ribose moiety
    uridine_complete_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](COP(O)(=O)O)O[C@H](n2ccc(=O)[nH]c2=O)[C@H]1O")
    
    # UDP linkage pattern (uridine + diphosphate)
    udp_linkage_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](COP([O-])(=O)OP(O)(=O)[O-])O[C@H](n2ccc(=O)[nH]c2=O)[C@H]1O")

    # Sugar moiety pattern can be defined using generic sugar scaffold
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O[C@@H]1*)")  # Flexible to cover sugars
    
    # Check for complete uridine pattern
    if not mol.HasSubstructMatch(uridine_complete_pattern):
        return False, "No complete uridine moiety"

    # Check for specific UDP linkage pattern
    if not mol.HasSubstructMatch(udp_linkage_pattern):
        return False, "No UDP linkage found"

    # Check for sugar component overlap with UDP linkage
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No recognizable sugar moiety linking to UDP"
    
    return True, "Contains specific uridine and diphosphate linkage typical of UDP-sugars with distinctive sugar moiety"