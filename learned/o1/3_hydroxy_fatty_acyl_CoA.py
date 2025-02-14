"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA results from the condensation of Coenzyme A (CoA) with a 3-hydroxy fatty acid via a thioester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the CoA moiety (simplified)
    coa_smarts = Chem.MolFromSmarts("""
        [#16X2]-[#6]-[#6]-[#7]-[#6](=O)-[#6]-[#7]-[#6](=O)-[#6H]-[#8]-[#6]([#6])([#6])-[#6]-[#8]-[P](=O)([O-])-[O]-[P](=O)([O-])-
        [O]-[#6]-[#8]-[#6H]-1-[#8]-[#6H]-[#6H]-1-[O]-[P](=O)([O-])[O-]-c1nc2c(n1)[nH]cnc2N
    """)
    
    if not mol.HasSubstructMatch(coa_smarts):
        return False, "CoA moiety not found"

    # Define a SMARTS pattern for the thioester linkage connecting CoA to acyl chain
    thioester_smarts = Chem.MolFromSmarts('C(=O)S[#16X2]-[#6]-[#6]-[#7]-[#6](=O)-[#6]')
    thioester_matches = mol.GetSubstructMatches(thioester_smarts)
    if not thioester_matches:
        return False, "No thioester linkage connecting CoA to acyl chain found"

    # Iterate over thioester matches to find the acyl chain
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon
        sulfur_idx = match[2]      # Sulfur atom (S in thioester linkage)

        # From carbonyl carbon, identify the acyl chain
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)

        # Identify the C2 carbon (attached to carbonyl carbon, not sulfur)
        c2_atom = None
        for neighbor in carbonyl_c.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != sulfur_idx:
                c2_atom = neighbor
                break
        if c2_atom is None:
            continue  # No acyl chain found

        # Identify C3 atom (next carbon in the chain)
        c3_atom = None
        for neighbor in c2_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != carbonyl_c_idx:
                c3_atom = neighbor
                break
        if c3_atom is None:
            continue  # Acyl chain too short

        # Check for hydroxy group on C3 atom
        has_hydroxy = False
        for neighbor in c3_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                # Ensure it's a hydroxy group (oxygen with only one bond)
                if neighbor.GetDegree() == 1:
                    has_hydroxy = True
                    break
        if has_hydroxy:
            return True, "Contains 3-hydroxy group on acyl chain attached to CoA via thioester linkage"

    return False, "Does not contain 3-hydroxy group on acyl chain attached to CoA via thioester linkage"