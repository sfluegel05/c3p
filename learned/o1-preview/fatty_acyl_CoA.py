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

    # Define Coenzyme A substructure pattern (simplified)
    # Focused on key functional groups: adenine ring, ribose sugar, diphosphate, and pantetheine unit
    coa_smarts = 'NC(=O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc32)[C@@H](O)[C@H]1O'

    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Error in Coenzyme A SMARTS pattern"

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # Define thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts('C(=O)SCCNC(=O)')
    if thioester_pattern is None:
        return False, "Error in thioester SMARTS pattern"

    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage to CoA found"

    # Identify the acyl chain attached via the thioester
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Get the carbonyl carbon of the thioester linkage
    # Assume first match
    match = thioester_matches[0]
    carbonyl_c_idx = match[0]  # Index of the carbonyl carbon

    # Traverse the acyl chain starting from the carbonyl carbon
    visited = set(match)  # Start with atoms in the thioester_pattern to avoid revisiting
    def traverse_acyl_chain(atom_idx):
        chain_atoms = []
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIdx() in visited:
            return chain_atoms
        visited.add(atom.GetIdx())
        if atom.GetAtomicNum() != 6:
            return chain_atoms  # Only consider carbon atoms
        chain_atoms.append(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in visited:
                continue
            bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
            if bond.GetBondType() in [Chem.BondType.SINGLE, Chem.BondType.DOUBLE]:
                chain_atoms.extend(traverse_acyl_chain(neighbor_idx))
        return chain_atoms

    acyl_chain_atoms = traverse_acyl_chain(carbonyl_c_idx)
    # Exclude the carbonyl carbon to get the length of the hydrocarbon chain
    acyl_chain_length = len(acyl_chain_atoms) - 1
    if acyl_chain_length < 4:
        return False, f"Acyl chain length is {acyl_chain_length}, too short for fatty acid"

    return True, f"Contains fatty acyl-CoA with acyl chain length {acyl_chain_length}"