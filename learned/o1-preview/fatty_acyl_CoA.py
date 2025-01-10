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

    # Define a more general Coenzyme A substructure pattern
    # Capture key features: adenine ring, ribose sugar, diphosphate, and pantetheine unit
    coa_smarts = '[#8]-[#3]-[#15](=O)([*])-[#8]-[#15](=O)([*])-[#8]-[#6]-1-[#8]-[#6](-[#8])-[#6](-[#8])-O1-[#7]-2-[#6]-3(=[#7]-[#6](=[#7]-[#6]-2=[#7]-3)-[#7])'

    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Error in Coenzyme A SMARTS pattern"

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # Define thioester linkage pattern connecting fatty acid to CoA
    thioester_pattern = Chem.MolFromSmarts('C(=O)SCCN')
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
    visited = set()
    def traverse_acyl_chain(atom_idx):
        chain_atoms = []
        stack = [atom_idx]
        while stack:
            idx = stack.pop()
            if idx in visited:
                continue
            visited.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # Only consider carbon atoms
            chain_atoms.append(idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in visited:
                    continue
                bond = mol.GetBondBetweenAtoms(idx, neighbor_idx)
                # Stop traversal if we reach the sulfur atom (to avoid including CoA carbons)
                if neighbor.GetAtomicNum() == 16:
                    continue
                # Exclude atoms towards CoA (nitrogen connected to the sulfur)
                if neighbor.GetAtomicNum() == 7 and neighbor_idx in match:
                    continue
                stack.append(neighbor_idx)
        return chain_atoms

    acyl_chain_atoms = traverse_acyl_chain(carbonyl_c_idx)
    # Exclude the carbonyl carbon to get the length of the hydrocarbon chain
    acyl_chain_length = len(acyl_chain_atoms) - 1
    if acyl_chain_length < 4:
        return False, f"Acyl chain length is {acyl_chain_length}, too short for fatty acid"

    return True, f"Contains fatty acyl-CoA with acyl chain length {acyl_chain_length}"