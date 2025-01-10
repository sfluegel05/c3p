"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: decanoate ester
Definition: A fatty acid ester resulting from the formal condensation of the carboxy group 
of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group pattern (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"

    # Pattern for decanoyl group with exact chain length and no branching
    # [CH3] ensures terminal methyl
    # Each carbon in chain specified as [CH2] to prevent branching
    # ($) ensures no additional connections beyond what's specified
    decanoyl_pattern = Chem.MolFromSmarts("""
        [OX2][CX3](=[OX1])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3]
    """)
    
    if not mol.HasSubstructMatch(decanoyl_pattern):
        return False, "No unbranched decanoyl group found"
    
    # Get matches
    matches = mol.GetSubstructMatches(decanoyl_pattern)
    
    # Check each potential match
    for match in matches:
        # Verify atoms are not in ring
        if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
            continue
            
        # Verify the ester oxygen is connected to carbon (not part of phosphate, etc)
        o_atom = mol.GetAtomWithIdx(match[0])
        neighbors = [x for x in o_atom.GetNeighbors() if x.GetIdx() != match[1]]
        if len(neighbors) != 1 or neighbors[0].GetAtomicNum() != 6:
            continue
            
        return True, "Contains decanoyl group connected via ester linkage with exactly 10 carbons"
        
    return False, "No valid decanoate ester group found"