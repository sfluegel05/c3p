"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA (greater than C22) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA substructure
    coa_substructure = Chem.MolFromSmiles('SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_substructure):
        return False, "CoA substructure not found"

    # Define the ester/thioester substructure (C(=O)S) to identify the carbonyl carbon connected to the fatty acid
    carbonyl_substructure = Chem.MolFromSmarts('C(=O)S')
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_substructure)

    if len(carbonyl_matches) == 0:
        return False, "No carbonyl linked to sulfur found"
    
    if len(carbonyl_matches) > 1:
         return False, "Multiple carbonyls linked to sulfur found, cannot identify fatty acid chain"

    # Get the carbon atom index in the carbonyl group
    carbonyl_carbon_index = carbonyl_matches[0][0]

    # trace the fatty acid chain by finding neighboring carbon atoms, excluding the carbonyl carbon
    fatty_acid_carbon_count = 0
    current_atom_index = carbonyl_carbon_index
    visited_atoms = {current_atom_index}

    while True:
        found_next_carbon = False
        for neighbor in mol.GetAtomWithIdx(current_atom_index).GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited_atoms:
               fatty_acid_carbon_count += 1
               visited_atoms.add(neighbor.GetIdx())
               current_atom_index = neighbor.GetIdx()
               found_next_carbon = True
               break

        if not found_next_carbon:
            break

    if fatty_acid_carbon_count > 22:
       return True, f"Fatty acid chain contains {fatty_acid_carbon_count} carbons, which is greater than 22"
    else:
      return False, f"Fatty acid chain contains {fatty_acid_carbon_count} carbons, which is not greater than 22"