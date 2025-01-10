"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: long-chain fatty alcohol
Definition: A fatty alcohol with a chain length ranging from C13 to C22
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13:
        return False, f"Insufficient carbons (C{c_count}, need C13-C22)"
    if c_count > 22:
        return False, f"Too many carbons (C{c_count}, maximum is C22)"

    # Look for alcohol groups (exclude phenols)
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    phenol_pattern = Chem.MolFromSmarts("c[OH]")
    
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No aliphatic alcohol group found"
    
    # Count alcohol groups
    alcohol_matches = len(mol.GetSubstructMatches(alcohol_pattern))
    if alcohol_matches > 3:  # Allow up to 3 alcohol groups (some examples have diols/triols)
        return False, f"Too many alcohol groups ({alcohol_matches})"

    # Exclude carboxylic acids as main feature
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxyl_pattern):
        # Check if there's still an alcohol group after excluding the acid
        acid_o_count = len(mol.GetSubstructMatches(carboxyl_pattern))
        if alcohol_matches <= acid_o_count:
            return False, "Primary feature is carboxylic acid, not alcohol"

    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 2:  # Allow up to 2 rings (some natural products have rings)
        return False, f"Too many rings ({ring_count})"

    # Check for excessive heteroatoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if n_count + s_count + p_count > 2:  # Allow up to 2 heteroatoms
        return False, "Too many heteroatoms for a fatty alcohol"

    # Calculate fraction of carbons and oxygens to ensure mainly hydrocarbon nature
    total_atoms = mol.GetNumAtoms()
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count / total_atoms < 0.6:  # At least 60% should be carbon
        return False, "Not primarily hydrocarbon in nature"
        
    return True, f"Contains alcohol group with appropriate carbon count (C{c_count})"