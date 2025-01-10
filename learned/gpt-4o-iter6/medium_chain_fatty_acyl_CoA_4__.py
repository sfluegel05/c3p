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

    # Look for characteristic CoA moiety
    coa_smarts = "[cR2]1ncnc2n(cnc12)[C@@H]3O[C@H](COP([O-])(=O)O)C[C@H](O)[C@H]3O"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA adenylate pattern detected"

    # Check for the thioester linkage C(=O)S part
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester bond found"

    # Verify medium-chain fatty acid length (evaluate linear chains starting from carbon after thioester)
    valid_fatty_acid = False
    for match in thioester_matches:
        thioester_carbon_idx = match[0]  # C bound to thioester
        neighbors = mol.GetAtomWithIdx(thioester_carbon_idx).GetNeighbors()
        
        # Find the carbon atom in the acyl chain
        acyl_chain_carbon = next((n for n in neighbors if n.GetAtomicNum() == 6), None)

        if acyl_chain_carbon:
            # Use breadth-first search to traverse linear acyl chain
            chain_length = 0
            visited = set()
            to_visit = [acyl_chain_carbon.GetIdx()]
            while to_visit:
                atom_idx = to_visit.pop(0)
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():  # Non-aromatic carbon
                    chain_length += 1
                    visited.add(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                            to_visit.append(neighbor.GetIdx())

            if 6 <= chain_length <= 12:
                valid_fatty_acid = True
                break

    if not valid_fatty_acid:
        return False, "Fatty acid chain not within medium-chain length (6 to 12 non-aromatic carbons)"

    # Confirm the deprotonated state to carry a -4 charge
    total_charge = rdmolops.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Total charge is not -4, found {total_charge}"

    return True, "Molecule matches medium-chain fatty acyl-CoA(4-) structure"