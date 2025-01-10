"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
"""
Classifies: aliphatic aldoxime
Definition: Any aldoxime derived from an aliphatic aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for aldoxime group pattern [H]C=NOH or C=NOH
    # Multiple SMARTS patterns to catch different representations
    aldoxime_patterns = [
        Chem.MolFromSmarts("[CH1](=N[OH1])"), # Explicit hydrogen
        Chem.MolFromSmarts("[CH0](=N[OH1])"),  # For cases where H is on other side
        Chem.MolFromSmarts("[CH1](=NO)"),      # Implicit H on oxygen
        Chem.MolFromSmarts("[CH0](=NO)")       # Alternative representation
    ]
    
    found_aldoxime = False
    for pattern in aldoxime_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_aldoxime = True
            break
            
    if not found_aldoxime:
        return False, "No aldoxime group (HC=NOH) found"

    # Check for aromatic atoms
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic atoms - must be aliphatic"
        
    # Get the carbon atom that's part of the C=N bond
    c_n_pattern = Chem.MolFromSmarts("[#6]=[#7][OH1]")
    if c_n_pattern is None:
        return False, "Error in SMARTS pattern"
        
    matches = mol.GetSubstructMatches(c_n_pattern)
    if not matches:
        return False, "No C=N-OH group found"
        
    # Check the hybridization of the carbon atom
    for match in matches:
        c_atom = mol.GetAtomWithIdx(match[0])
        if c_atom.GetHybridization() != Chem.HybridizationType.SP2:
            return False, "Carbon in C=N must be sp2 hybridized"
            
        # Check that the carbon has exactly one hydrogen
        # (characteristic of aldoxime vs ketoxime)
        if c_atom.GetTotalNumHs() != 1:
            return False, "Carbon must have exactly one hydrogen (aldoxime)"
            
    return True, "Contains aliphatic aldoxime group (HC=NOH)"