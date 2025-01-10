"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
A sterol ester obtained by formal condensation of the carboxy group of any 
carboxylic acid with the 3-hydroxy group of cholesterol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define multiple SMARTS patterns for steroid core to catch different representations
    steroid_patterns = [
        # Basic steroid core with flexible bond types
        "[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6](~[#6]~2)~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~3~1",
        # Alternative pattern focusing on ring connectivity
        "C1CC2CCC3(C)C(CCC4C3CC=C4)C2CC1",
        # More specific pattern including double bond
        "C1CC2=CCC3C(CC[C@@H]4CCCC[C@]34C)C2(C)CC1"
    ]
    
    core_found = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            core_found = True
            break
            
    if not core_found:
        return False, "No steroid core structure found"

    # Look for ester group at position 3 (-O-C(=O)-)
    # More specific pattern that includes connection to ring
    ester_pattern = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6;X4]~[#6]~[#6]~[#6]~1[OX2][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Missing ester group in correct position"

    # Check for characteristic cholesterol features
    
    # Double bond (typically at C5-C6, but allow flexibility)
    double_bond_patterns = [
        "C=CC1CCC2C(C1)CCC",
        "CC1CC=C2CC1CCC2"
    ]
    double_bond_found = False
    for pattern in double_bond_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            double_bond_found = True
            break
    
    if not double_bond_found:
        return False, "Missing characteristic steroid double bond"

    # Check for branched side chain with flexible pattern
    side_chain_patterns = [
        "CCC(C)CCC(C)C",
        "[#6][#6][#6]([#6])[#6][#6][#6]([#6])[#6]"
    ]
    side_chain_found = False
    for pattern in side_chain_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            side_chain_found = True
            break
            
    if not side_chain_found:
        return False, "Missing characteristic cholesterol side chain"

    # Verify oxygen count (exactly 2 for ester group)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count != 2:
        return False, f"Found {oxygen_count} oxygens, should be exactly 2 for cholesteryl ester"

    # Count carbons - cholesteryl esters typically have >27 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27:
        return False, "Too few carbons for cholesteryl ester"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:  # Adjusted threshold
        return False, "Molecular weight too low for cholesteryl ester"

    return True, "Contains cholesterol core with ester-linked fatty acid chain"