"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: wax ester
A fatty acid ester resulting from the condensation of the carboxy group of a fatty acid 
with the alcoholic hydroxy group of a fatty alcohol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    
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

    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Count carbons, oxygens, and check for other elements
    c_count = 0
    o_count = 0
    other_elements = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            c_count += 1
        elif atom.GetAtomicNum() == 8:  # Oxygen
            o_count += 1
        elif atom.GetAtomicNum() not in [1, 6, 8]:  # Not H, C, or O
            other_elements = True
            
    if other_elements:
        return False, "Contains elements other than C, H, and O"
    
    if o_count != 2:
        return False, "Must have exactly 2 oxygens (one ester group)"
        
    if c_count < 20:
        return False, "Carbon chain too short for wax ester (need C20+ total)"

    # Check for long carbon chains on both sides of ester
    # First split molecule at ester
    fragments = Chem.FragmentOnBonds(mol, [bond.GetIdx() for bond in mol.GetBonds() 
                                         if bond.GetBeginAtom().GetAtomicNum() == 8 
                                         and bond.GetEndAtom().GetAtomicNum() == 6])
    if fragments is None:
        return False, "Could not analyze chain lengths"
        
    # Count carbons in each fragment
    fragment_sizes = []
    for frag in Chem.GetMolFrags(fragments, asMols=True):
        c_in_frag = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        fragment_sizes.append(c_in_frag)
        
    if len(fragment_sizes) < 2:
        return False, "Could not identify two distinct chains"
        
    min_chain = min(fragment_sizes)
    if min_chain < 8:
        return False, f"Shortest chain has only {min_chain} carbons, need 8+"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Too few rotatable bonds for wax ester"

    # Look for problematic groups that would make it not a wax ester
    problematic_groups = [
        "[OH]",  # alcohol (except the one in ester)
        "[CX3](=O)[OX2H]",  # carboxylic acid
        "[NX3]",  # amine
        "[SX2]",  # thiol/sulfide
    ]
    
    for smart in problematic_groups:
        pattern = Chem.MolFromSmarts(smart)
        if mol.HasSubstructMatch(pattern):
            return False, "Contains additional functional groups"

    return True, "Contains one ester group connecting two long hydrocarbon chains"