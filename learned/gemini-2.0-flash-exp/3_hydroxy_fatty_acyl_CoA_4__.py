"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA(4-) consists of a fatty acid chain with a 3-hydroxy group,
    linked to coenzyme A via a thioester, with 4 negative charges on the phosphate groups

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core CoA Pattern (without charges)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA core structure found."

    # Check for the thioester link
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester link found."

    # Check for the 3-hydroxy group (relaxed pattern)
    hydroxy_pattern = Chem.MolFromSmarts("[CX4][CX4](O)")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)

    # Check if a 3-hydroxy group is within a fatty acyl chain
    found_valid_3_hydroxy = False
    for match in hydroxy_matches:
        # Get the index of the carbon with OH
        oh_carbon_index = match[1]
        oh_carbon = mol.GetAtomWithIdx(oh_carbon_index)
        
        #Check neighbors
        neighbors = [n for n in oh_carbon.GetNeighbors() if n.GetAtomicNum() == 6] # Only carbon neighbors
        if len(neighbors) != 2:
            continue

        # Check if neighbors are connected to chain of carbons
        c1, c2 = neighbors
        c1_neighbors = [n for n in c1.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != oh_carbon_index]
        c2_neighbors = [n for n in c2.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != oh_carbon_index]

        if len(c1_neighbors) > 0 and len(c2_neighbors) > 0:
              #check that the number of carbons is at least 3, not counting the hydroxyl carbon
              chain_carbons = 0
              visited_atoms = {oh_carbon_index}
              
              def count_chain_recursive(atom, depth):
                  nonlocal chain_carbons
                  if depth > 5:
                      return
                  
                  visited_atoms.add(atom.GetIdx())
                  for neighbor in atom.GetNeighbors():
                      if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited_atoms:
                           chain_carbons += 1
                           count_chain_recursive(neighbor, depth+1)

              count_chain_recursive(c1, 0)
              count_chain_recursive(c2, 0)

              if chain_carbons >= 3:
                   found_valid_3_hydroxy = True
                   break

    if not found_valid_3_hydroxy:
        return False, "No 3-hydroxy group found or not on fatty acid chain."
    
    # Check for four negative charges on phosphates.
    charge = 0
    for atom in mol.GetAtoms():
        charge += atom.GetFormalCharge()
    if charge != -4:
        return False, f"Incorrect overall charge ({charge}), must be -4."
    

    # Additional checks:
    # Counting the specific atoms - check they are within the expected ranges,
    # and that there is at least 1 sulfur and 4 phosphorus.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if s_count < 1:
        return False, "Must have at least 1 sulfur"
    if p_count != 4:
        return False, "Must have 4 phosphorus atoms"
    if c_count < 20:
       return False, "Too few carbon atoms to be a fatty acyl CoA"

    return True, "Molecule matches the criteria for a 3-hydroxy fatty acyl-CoA(4-)."