"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA

"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA is an acyl-CoA resulting from the condensation of the thiol group 
    of coenzyme A with the carboxy group of a fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define Coenzyme A substructure pattern
    # The pattern includes adenine, ribose, diphosphate, pantetheine moiety
    coa_smarts = """
    O=P(O)(O)OP(O)(O)OC[C@H]1O[C@H]([C@@H](O)[C@H]1O)n2cnc3c(ncnc32)N
    """

    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Error in Coenzyme A SMARTS pattern"

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # Define thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage to CoA found"

    # Now identify the acyl chain attached via the thioester
    # Find the carbonyl carbon of the thioester
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Get the acyl chain starting from the carbonyl carbon
    # Assume first match
    match = thioester_matches[0]
    carbonyl_c_idx = match[0]  # Index of the carbonyl carbon

    # Traverse the acyl chain
    visited = set()
    def traverse_acyl_chain(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIdx() in visited:
            return 0
        visited.add(atom.GetIdx())
        if atom.GetAtomicNum() != 6:
            return 0  # Only count carbon atoms
        length = 1
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
            if bond.GetBondType() not in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
                continue  # Ignore non-single/double bonds
            if neighbor_idx in match:
                continue  # Avoid going back to known parts of the thioester
            length += traverse_acyl_chain(neighbor_idx)
        return length

    acyl_length = traverse_acyl_chain(carbonyl_c_idx) - 1  # Subtract 1 to exclude carbonyl carbon
    if acyl_length < 4:
        return False, f"Acyl chain length is {acyl_length}, too short for fatty acid"

    return True, f"Contains fatty acyl-CoA with acyl chain length {acyl_length}"