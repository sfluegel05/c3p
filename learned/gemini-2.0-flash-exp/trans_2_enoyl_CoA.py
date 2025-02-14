"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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
    
    # Define the thioester group
    thioester_pattern = Chem.MolFromSmarts("S[CX3](=O)")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
         return False, f"Found {len(thioester_matches)} thioester groups, need exactly 1"
    
    # Define the trans double bond pattern
    trans_double_bond_pattern = Chem.MolFromSmarts("[CX4]=[CX4]")
    trans_double_bond_matches = mol.GetSubstructMatches(trans_double_bond_pattern)
    
    # Check for trans configuration.
    found_correct_bond = False
    for match in trans_double_bond_matches:
        bond = mol.GetBondBetweenAtoms(match[0],match[1])
        if bond.GetStereo() != Chem.BondStereo.STEREOE:
           continue # skip if not trans
        
        #Find atoms connected to the double bond
        atom1 = mol.GetAtomWithIdx(match[0])
        atom2 = mol.GetAtomWithIdx(match[1])

        # find the carbonyl carbon
        carbonyl_carbon_idx = None
        for neighbor in atom1.GetNeighbors():
            if neighbor.GetIdx() not in match:
              neighbor_atom = mol.GetAtomWithIdx(neighbor.GetIdx())
              if neighbor_atom.GetAtomicNum() == 8 and neighbor_atom.GetTotalValence() == 2:
                   # This is the carbonyl oxygen, so skip
                   continue
              elif neighbor_atom.GetAtomicNum() == 8 and neighbor_atom.GetTotalValence() == 1:
                   # not a carbonyl
                   continue
              elif neighbor_atom.GetAtomicNum() == 6 and neighbor_atom.GetTotalValence() == 3:
                  # this is the thioester carbonyl carbon
                  carbonyl_carbon_idx = neighbor.GetIdx()
        
        if carbonyl_carbon_idx is None:
          continue # skip if the double bond does not come from an acyl chain
        
        # find the thioester sulfur
        thio_sulfur_idx = None
        for neighbor in mol.GetAtomWithIdx(carbonyl_carbon_idx).GetNeighbors():
            neighbor_atom = mol.GetAtomWithIdx(neighbor.GetIdx())
            if neighbor_atom.GetAtomicNum() == 16:
                 thio_sulfur_idx = neighbor.GetIdx()
                 break
        
        if thio_sulfur_idx is None:
           continue # skip if no sulfur connected

        # find if the sulfur belongs to the CoA
        for sulfur_match in mol.GetSubstructMatches(Chem.MolFromSmarts('S')):
            if sulfur_match[0] == thio_sulfur_idx:
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