"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the CoA adenylate part and general structure indicative of CoA
    adenylate_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2N")
    coA_general_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)O")
    if not (mol.HasSubstructMatch(adenylate_pattern) and mol.HasSubstructMatch(coA_general_pattern)):
        return False, "No complete Coenzyme A structure found"

    # Check for the thioester linkage C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) == 0:
        return False, "No thioester bond found"

    # Verify medium-chain fatty acid length (evaluate linear chains starting from carbon after thioester)
    valid_fatty_acid = False
    for match in thioester_matches:
        thioester_carbon_idx = match[1]  # S bound to thioester carbon
        neighbors = mol.GetAtomWithIdx(thioester_carbon_idx).GetNeighbors()
        
        # Find the carbon atom in the acyl chain
        next_c_atom = None
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                next_c_atom = neighbor.GetIdx()
                break

        if next_c_atom is not None:
            chain_length = 0
            visited = set()

            # Perform Breadth-First Search (BFS) to traverse carbon chain
            to_visit = [next_c_atom]
            while to_visit:
                atom_idx = to_visit.pop(0)
                if atom_idx in visited:
                    continue
                if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6:
                    chain_length += 1
                visited.add(atom_idx)
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                        to_visit.append(neighbor.GetIdx())

            if 6 <= chain_length <= 12:
                valid_fatty_acid = True
                break

    if not valid_fatty_acid:
        return False, "Fatty acid chain not within medium-chain length (6 to 12 carbons)"

    # Confirm the deprotonated state to carry a -4 charge
    total_charge = rdmolops.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Total charge is not -4, found {total_charge}"

    return True, "Molecule matches medium-chain fatty acyl-CoA(4-) structure"