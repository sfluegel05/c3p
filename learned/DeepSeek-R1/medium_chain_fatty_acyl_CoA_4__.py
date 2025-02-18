"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:XXXXX medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if matches criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Verify CoA core structure components
    # Match adenine-ribose-phosphate pattern
    coa_core = Chem.MolFromSmarts(
        "N1C=NC2=C1N=CN2[C@H]1[C@@H]([C@@H]([C@H](O1)COP(OP(OC[C@H]2O[C@H]([C@H]([C@@H]2O)O)CO)([O-])=O)([O-])=O)([O-])=O)O"
    )
    if not mol.HasSubstructMatch(coa_core):
        return False, "Missing CoA core structure"

    # Check for thioester-linked acyl group (S-C(=O)-R)
    thioester = Chem.MolFromSmarts("[SX2][CX3](=O)")
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester bond found"

    # Get the acyl chain length
    try:
        sulfur = thioester_matches[0][0]
        carbonyl = thioester_matches[0][1]
    except IndexError:
        return False, "Invalid thioester bond"

    # Traverse the acyl chain
    chain_atoms = set()
    stack = [carbonyl]
    while stack:
        atom = stack.pop()
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() == sulfur:
                continue
            if neighbor.GetAtomicNum() == 1: