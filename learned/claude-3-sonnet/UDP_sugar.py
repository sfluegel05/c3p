"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: UDP-sugar
A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached 
to an unspecified sugar via an anomeric diphosphate linkage.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
        
    # Check for uracil ring - more flexible pattern
    uracil_pattern = Chem.MolFromSmarts("[#7]1[#6](=[O,S])[#7][#6][#6]1=[O,S]")
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "No uracil ring found"
    
    # Check for ribose attached to uracil (uridine) - more flexible pattern
    uridine_pattern = Chem.MolFromSmarts("[#7]1[#6](=[O,S])[#7][#6][#6]1=[O,S].[#6]1[#8][#6][#6][#6][#6]1")
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "No uridine moiety found"
        
    # Check for diphosphate group - more flexible pattern
    diphosphate_pattern = Chem.MolFromSmarts("[O,S]-[P](=[O,S])([O,S])-[O,S]-[P](=[O,S])([O,S])[O,S]")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate linkage found"
        
    # Check for sugar characteristics (ring with multiple oxygens)
    # More flexible pattern that can match various sugar modifications
    sugar_pattern = Chem.MolFromSmarts("[#6]1[#8,#7][#6][#6][#6][#6]1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"
        
    # Verify UDP core structure
    udp_core = Chem.MolFromSmarts("""
        [#7]1[#6](=[O,S])[#7][#6][#6]1=[O,S]  # Uracil
        .[#6]1[#8][#6][#6][#6][#6]1            # Ribose
        .[O,S]-[P](=[O,S])([O,S])-[O,S]-[P]    # Diphosphate
    """)
    if not mol.HasSubstructMatch(udp_core):
        return False, "Missing UDP core structure"
    
    # Count phosphorus atoms - should be exactly 2 for UDP
    p_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[P]")))
    if p_count != 2:
        return False, f"Found {p_count} phosphorus atoms, need exactly 2"
    
    # Check connectivity between components
    # Look for ribose-phosphate and phosphate-sugar connections
    connections_pattern = Chem.MolFromSmarts("""
        [#6]1[#8][#6][#6][#6][#6]1[#8][P]  # Ribose-phosphate
        .[P][#8][#6]                        # Phosphate-sugar
    """)
    if not mol.HasSubstructMatch(connections_pattern):
        return False, "Components not properly connected"
        
    return True, "Contains UDP moiety connected to sugar via anomeric diphosphate linkage"