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
    A fatty acyl-CoA results from the condensation of the thiol group of coenzyme A
    with the carboxy group of any fatty acid.

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

    # Check for adenosine moiety
    adenosine_pattern = Chem.MolFromSmarts("n1c[nH]c2c1ncnc2")
    if not mol.HasSubstructMatch(adenosine_pattern):
        return False, "No adenosine moiety found"

    # Check for diphosphate bridge
    diphosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)O")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate bridge found"

    # Check for pantetheine moiety with thioester linkage
    pantetheine_thioester_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)C(O)C(C)(C)COP")
    if not mol.HasSubstructMatch(pantetheine_thioester_pattern):
        return False, "No pantetheine moiety found"

    # Check for fatty acyl chain attached via thioester linkage
    # Find the carbonyl carbon of the thioester
    carbonyl_carbons = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)S"))
    if not carbonyl_carbons:
        return False, "No thioester carbonyl carbon found"
    carbonyl_carbon_idx = carbonyl_carbons[0][0]

    # Traverse the fatty acyl chain
    fatty_chain_length = 0
    visited_atoms = set()
    atom_stack = [mol.GetAtomWithIdx(carbonyl_carbon_idx)]

    while atom_stack:
        atom = atom_stack.pop()
        atom_idx = atom.GetIdx()
        if atom_idx in visited_atoms:
            continue
        visited_atoms.add(atom_idx)

        if atom.GetAtomicNum() == 6:  # Carbon atom
            fatty_chain_length += 1

            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
                # Avoid going towards the sulfur atom of thioester linkage
                if neighbor.GetAtomicNum() == 16:
                    continue
                if neighbor_idx not in visited_atoms:
                    atom_stack.append(neighbor)

    if fatty_chain_length < 4:
        return False, f"Fatty acyl chain is too short ({fatty_chain_length} carbons)"

    return True, "Molecule is a fatty acyl-CoA with appropriate coenzyme A structure"