"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for adenine using multiple possible patterns
    adenine_patterns = [
        "c1ncnc2[nH]cnc12",  # Basic adenine core
        "c1ncnc2ncnc12",     # Alternative representation
        "c1nc(N)nc2[nH]cnc12",  # With amino group
        "c1nc(N)nc2ncnc12",     # Another amino form
        "[nH]1cnc2c(ncnc2n1)",  # Different tautomer
        "n1cnc2c(N)ncnc12"      # Yet another form
    ]
    
    found_adenine = False
    for pattern in adenine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_adenine = True
            break
    
    if not found_adenine:
        return False, "No adenine moiety found"

    # Look for phosphate groups characteristic of CoA
    phosphate_patterns = [
        "OP(O)(=O)OP(O)(=O)O",
        "P(O)(O)(=O)OP(O)(O)=O"
    ]
    
    found_phosphates = False
    for pattern in phosphate_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_phosphates = True
            break
            
    if not found_phosphates:
        return False, "Missing characteristic phosphate groups of CoA"
    
    # Check for pantetheine arm with thioester
    pantetheine_patterns = [
        "NCCC(=O)NCCS",
        "SCCNC(=O)CCNC"
    ]
    
    found_pantetheine = False
    for pattern in pantetheine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_pantetheine = True
            break
            
    if not found_pantetheine:
        return False, "Missing characteristic pantetheine arm of CoA"

    # Check for thioester linkage
    thioester_patterns = [
        "C(=O)S",
        "SC(=O)"
    ]
    
    found_thioester = False
    for pattern in thioester_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_thioester = True
            break
            
    if not found_thioester:
        return False, "No thioester linkage found"

    # Count number of double bonds in the fatty acid portion
    # Look for C=C bonds that are not part of the adenine ring
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    
    # Get aromatic atoms
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAromaticAtoms()}
    
    # Filter out double bonds where either carbon is aromatic
    non_aromatic_double_bonds = [
        bond for bond in double_bonds 
        if not (bond[0] in aromatic_atoms or bond[1] in aromatic_atoms)
    ]
    
    if len(non_aromatic_double_bonds) == 0:
        return False, "No carbon-carbon double bonds found in fatty acid portion"

    # Check for reasonable chain length (at least 4 carbons in fatty acid portion)
    carbon_chain_patterns = [
        "CCCC",
        "C~C~C~C"  # More flexible pattern allowing any bonds between carbons
    ]
    
    found_chain = False
    for pattern in carbon_chain_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_chain = True
            break
            
    if not found_chain:
        return False, "Fatty acid portion too short"

    # Success case - molecule has all required features
    return True, f"Contains CoA thioester linkage and {len(non_aromatic_double_bonds)} C=C double bonds in fatty acid portion"