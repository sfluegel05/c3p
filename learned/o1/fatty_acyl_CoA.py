"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA is an acyl-CoA that results from the formal condensation 
    of the thiol group of coenzyme A with the carboxy group of any fatty acid.

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

    # Check for thioester linkage (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Simplify adenosine moiety pattern (remove stereochemistry)
    adenosine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2")
    if not mol.HasSubstructMatch(adenosine_pattern):
        return False, "No adenosine moiety found, not coenzyme A"

    # Simplify pantetheine moiety pattern
    pantetheine_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)C(O)C(C)(C)COP")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine moiety found, not coenzyme A"

    # Check for full coenzyme A structure
    coenzymeA_pattern = Chem.MolFromSmarts("O=C(NCCSC=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@H]1O)n2cnc3c(n2)nc(N)[nH]c3=O")
    if not mol.HasSubstructMatch(coenzymeA_pattern):
        return False, "Full coenzyme A structure not found"

    # Check for fatty acyl chain attached via thioester linkage
    # The fatty acyl chain is connected to the carbonyl carbon of the thioester
    fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)[CH2][CH2]")
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)
    if not fatty_acyl_matches:
        return False, "No fatty acyl chain attached via thioester linkage found"

    # Determine length of fatty acyl chain (should be >4 carbons)
    # Start from the thioester carbonyl carbon and traverse the chain
    carbonyl_carbons = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)S"))
    if not carbonyl_carbons:
        return False, "No thioester carbonyl carbon found"
    carbonyl_carbon_idx = carbonyl_carbons[0][0]
    chain_length = 0
    visited = set()
    to_visit = [mol.GetAtomWithIdx(carbonyl_carbon_idx)]
    while to_visit:
        atom = to_visit.pop()
        if atom.GetIdx() in visited:
            continue
        visited.add(atom.GetIdx())
        if atom.GetAtomicNum() == 6:  # Carbon
            chain_length += 1
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    # Avoid traversing back to the coenzyme A moiety
                    if neighbor.HasSubstructMatch(Chem.MolFromSmarts("SCCNC(=O)")):
                        continue
                    to_visit.append(neighbor)
    if chain_length <= 4:
        return False, f"Fatty acyl chain is too short ({chain_length} carbons)"

    return True, "Contains coenzyme A moiety with fatty acyl chain attached via thioester linkage"