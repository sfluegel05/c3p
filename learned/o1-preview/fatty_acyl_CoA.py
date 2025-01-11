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
    A fatty acyl-CoA is an acyl-CoA that results from the formal condensation of the thiol group 
    of coenzyme A with the carboxy group of any fatty acid.

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
    
    # Define acyl-CoA linkage pattern: [NX3][CH2][CH2][SX2][CX3](=O)[C]
    acyl_coa_pattern = Chem.MolFromSmarts("[NX3][CH2][CH2][SX2][CX3](=O)[C]")
    if not mol.HasSubstructMatch(acyl_coa_pattern):
        return False, "No acyl-CoA linkage found"

    # Define adenine nucleotide moiety (adenine base)
    adenine_pattern = Chem.MolFromSmarts("n1c2ncnc2nc1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine moiety found"
    
    # Extract acyl chain length
    matches = mol.GetSubstructMatches(acyl_coa_pattern)
    if len(matches) == 0:
        return False, "No acyl-CoA linkage found"

    # Assuming the first match
    match = matches[0]
    acyl_carbon_idx = match[5]

    # Traverse the acyl chain starting from the acyl carbon
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
            if neighbor_idx == match[4]:  # Avoid going back to carbonyl carbon
                continue
            length += traverse_acyl_chain(neighbor_idx)
        return length

    acyl_length = traverse_acyl_chain(acyl_carbon_idx)

    if acyl_length < 4:
        return False, f"Acyl chain length is {acyl_length}, too short for fatty acid"

    return True, f"Contains acyl-CoA linkage with acyl chain length {acyl_length}"