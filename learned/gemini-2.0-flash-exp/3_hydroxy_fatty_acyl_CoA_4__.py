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
    # Trying to include the negative charges in the SMARTS query seems to fail
    # The charges will be handled later with a global analysis on the molecule
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA core structure found."

    # Check for the thioester link
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester link found."

    # Check for the 3-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[CX4][CX4](O)[CX4]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)

    # Extract positions of matches. Then check if one of them is within a fatty acyl chain
    found_valid_3_hydroxy = False
    for match in hydroxy_matches:
        # Get the index of the carbon with OH
        oh_carbon_index = match[1]
        oh_carbon = mol.GetAtomWithIdx(oh_carbon_index)
        if not oh_carbon.HasProp('_is_alkyl'):
            continue # skip if OH is not on an alkyl chain.
        # Check neighbours, at least two carbon atoms next to the OH carbon
        neighbours = [n.GetSymbol() for n in oh_carbon.GetNeighbors()]
        c_count = 0
        for n in oh_carbon.GetNeighbors():
            if n.GetAtomicNum() == 6:
                c_count += 1
                carbon_neighbours = [n2.GetSymbol() for n2 in n.GetNeighbors()]
                found_carbon_in_chain = False
                for nn in n.GetNeighbors():
                    if (nn.GetIdx() != oh_carbon.GetIdx()) and nn.GetAtomicNum()==6:
                        found_carbon_in_chain = True
                        break
                if not found_carbon_in_chain:
                    c_count = 0
                    break
                
        if c_count >= 2:
            found_valid_3_hydroxy = True
            break #found a 3-hydroxy group on the fatty acid chain
            
    if not found_valid_3_hydroxy:
        return False, "No 3-hydroxy group found or not on fatty acid chain."

    # Check for the presence of fatty acid chain
    # This is a loose check, so if fails there might still be a valid chain
    # The presence of a chain is enforced by the 3-hydroxy group
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
            return False, "Missing a fatty acid chain"
    
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