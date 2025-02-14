"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for CoA moiety
    coa_smarts = "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    
    if not mol.HasSubstructMatch(coa_pattern):
         return False, "Molecule does not contain the CoA moiety"

    # Find the sulfur atom in the CoA moiety
    match = mol.GetSubstructMatch(coa_pattern)
    sulfur_index = -1
    for atom_index in match:
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetAtomicNum() == 16:  # Sulfur is atomic number 16
            sulfur_index = atom_index
            break
    if sulfur_index == -1:
        return False, "Could not locate sulfur in CoA moiety"

    # Find the carbonyl carbon attached to the sulfur
    carbonyl_carbon_index = -1
    sulfur_atom = mol.GetAtomWithIdx(sulfur_index)
    for neighbor in sulfur_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and  neighbor.GetDoubleBondNeighbor() is not None and neighbor.GetDoubleBondNeighbor().GetAtomicNum() == 8: # Carbonyl carbon
            carbonyl_carbon_index = neighbor.GetIdx()
            break
    if carbonyl_carbon_index == -1:
        return False, "Could not locate carbonyl carbon next to sulfur"
    
    # Trace the fatty acid chain
    fatty_acid_carbon_chain = []
    current_atom = mol.GetAtomWithIdx(carbonyl_carbon_index)
    prev_atom = None

    while current_atom.GetAtomicNum() == 6: # keep going until the carbon chain ends
        fatty_acid_carbon_chain.append(current_atom.GetIdx())
        neighbors = current_atom.GetNeighbors()
        next_atom = None
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != (prev_atom.GetIdx() if prev_atom else -1):
              next_atom = neighbor
              break
        if next_atom:
           prev_atom = current_atom
           current_atom = next_atom
        else:
           break
    
    if len(fatty_acid_carbon_chain) < 12:
        return False, "Fatty acid chain too short"
        
    # Check saturation between 11th and 12th carbons
    carbon_11_index = fatty_acid_carbon_chain[10] if len(fatty_acid_carbon_chain) > 10 else -1
    carbon_12_index = fatty_acid_carbon_chain[11] if len(fatty_acid_carbon_chain) > 11 else -1
    
    if carbon_11_index == -1 or carbon_12_index == -1:
        return False, "Could not find carbons at positions 11 and 12"

    bond = mol.GetBondBetweenAtoms(carbon_11_index, carbon_12_index)
    if bond is None:
       return False, "No bond between carbon 11 and 12"
    
    if bond.GetBondType() != Chem.BondType.SINGLE:
        return False, "Bond between carbons 11 and 12 is not single (saturated)"
    
    return True, "The bond between carbons 11 and 12 is saturated in the fatty acyl chain"