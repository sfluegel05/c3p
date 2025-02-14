"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA is a coenzyme A molecule with a fatty acid chain
    attached via a thioester bond with a trans double bond at the 2-position
    relative to the carbonyl.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise.
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the CoA SMARTS pattern. This is the minimal core of the CoA moiety
    # we are checking to avoid issues with tautomers.
    coa_pattern = Chem.MolFromSmarts('SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
   
    # Define the trans-2-enoyl-CoA fragment with trans double bond
    # This SMARTS pattern specifically looks for a C=C with a trans config and with the acyl chain
    # and sulfur from coA attached.
    trans_2_enoyl_pattern = Chem.MolFromSmarts("[CX4][C@H]=[C@H][CX3](=O)[SX2]")
    
    matches = mol.GetSubstructMatches(trans_2_enoyl_pattern)

    found_correct_bond = False
    for match in matches:
      # check if the sulfur belongs to the coA using substruct match
      sulfur_idx = match[4]
      for sulfur_match in mol.GetSubstructMatches(Chem.MolFromSmarts('S')):
           if sulfur_match[0] == sulfur_idx:
               found_correct_bond = True
               break
      if found_correct_bond:
          break

    if not found_correct_bond:
        return False, "Could not find correct trans double bond"

    # Check molecular weight - generally > 700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for trans-2-enoyl-CoA"

    return True, "Molecule is a trans-2-enoyl-CoA"