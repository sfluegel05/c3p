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

    # More flexible CoA core pattern (ignores stereochemistry and exact connectivity)
    coa_core = Chem.MolFromSmarts(
        "[NH]C(=O)CCNC(=O)C(O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OC(C(O)C1O)n1cnc2c(N)ncnc12"
    )
    if not mol.HasSubstructMatch(coa_core):
        return False, "Missing CoA core structure"

    # Find thioester group (S-C(=O)-R)
    thioester = Chem.MolFromSmarts("[SX2][CX3](=O)")
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester bond found"

    # Get the sulfur and carbonyl atoms from the first thioester match
    sulfur_idx = thioester_matches[0][0]
    carbonyl_idx = thioester_matches[0][1]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)

    # Function to calculate longest carbon chain from starting atom
    def find_longest_chain(atom, visited=None):
        if visited is None:
            visited = set()
        if atom.GetIdx() in visited or atom.GetAtomicNum() != 6:
            return 0
        visited.add(atom.GetIdx())
        max_length = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() == sulfur_idx:  # Skip back to CoA
                continue
            # Skip oxygen/nitrogen neighbors (stay in hydrocarbon chain)
            if neighbor.GetAtomicNum() not in {6, 1}:
                continue
            # Avoid double-counting with bond type check
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetAtomicNum() != 6:
                continue
            length = find_longest_chain(neighbor, visited.copy())
            if length > max_length:
                max_length = length
        return 1 + max_length

    # Calculate chain length starting from carbonyl carbon (subtract 1 to exclude carbonyl itself)
    chain_length = find_longest_chain(carbonyl_atom) - 1
    if not (6 <= chain_length <= 12):
        return False, f"Main acyl chain length {chain_length} not in medium range (6-12)"

    # Check total molecular charge is -4
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != -4:
        return False, f"Total charge {total_charge} ≠ -4"

    return True, "Medium-chain fatty acyl-CoA(4-) with correct structure and charge"