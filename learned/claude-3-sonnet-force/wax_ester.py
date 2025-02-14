"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: CHEBI:38219 wax ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of the carboxy group
    of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Look for fatty acid chain (long carbon chain attached to carbonyl carbon)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing fatty acid chain"
    
    # Look for fatty alcohol chain (long carbon chain with terminal oxygen attached to ester oxygen)
    fatty_alcohol_pattern = Chem.MolFromSmarts("[OX2][CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_alcohol_matches = mol.GetSubstructMatches(fatty_alcohol_pattern)
    if not fatty_alcohol_matches:
        return False, "Missing fatty alcohol chain"

    # Check for minimum rotatable bonds (10) to ensure long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty chains"

    # Check for minimum carbon atoms (20)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for wax ester"

    # Check for exactly two oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 2:
        return False, "Must have exactly 2 oxygens (1 ester group)"

    return True, "Contains fatty acid ester with a fatty acid chain and a fatty alcohol chain"