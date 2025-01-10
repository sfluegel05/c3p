"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of the carboxy group of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ester functional group pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CX4H2][CX4H2]*")
    
    # Check if the molecule contains an ester group
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, "No ester group found or structure too complex"

    # Check the length of carbon chains flanking the ester group
    ester_idx = ester_matches[0]
    
    # Define the pattern for a long aliphatic chain (at least 10 carbons)
    long_chain_pattern = Chem.MolFromSmarts("[CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2]*")
    
    # Split the molecule to check chains on both sides of the ester
    ester_bond = [(ester_idx[1], ester_idx[2])] # potential ester bond
    Chem.Kekulize(mol, True)
    
    # Check for a significant long chain on each side
    found_left_chain, found_right_chain = False, False
    
    for chain in mol.GetSubstructMatches(long_chain_pattern):
        if ester_idx[0] in chain:
            found_left_chain = True
        if ester_idx[2] in chain:
            found_right_chain = True
    
    if not (found_left_chain and found_right_chain):
        return False, "One or both chains are too short"

    # The molecule should not be cyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings, unlikely to be a simple wax ester"

    # Check for unsaturation
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern)
    
    if has_double_bond:
        reason = "Contains single ester group with long carbon chains including double bonds - likely an unsaturated wax ester"
    else:
        reason = "Contains single ester group with long carbon chains - likely a saturated wax ester"
    
    return True, reason