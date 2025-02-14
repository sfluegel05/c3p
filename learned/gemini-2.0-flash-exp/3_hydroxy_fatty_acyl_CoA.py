"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA has a coenzyme A moiety linked via a thioester to the carbonyl of a 3-hydroxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    #  CoA pattern including the terminal phosphate
    coa_pattern = Chem.MolFromSmarts("COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A substructure not found."


    # Combined pattern for thioester and 3-hydroxy carbon (C=O-S-C-C(O))
    hydroxy_thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])S[CX4][CX4]([OX2])")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_thioester_pattern)
    if len(hydroxy_matches) == 0:
       return False, "No thioester with a 3-hydroxy group found"

    found_chain = False
    for match in hydroxy_matches:
       
        thioester_carbonyl_idx = match[0]  # C=O of thioester
        alpha_carbon_idx = match[2]      # carbon attached to S
        hydroxy_carbon_idx = match[3]  # The carbon with the OH group

        # Check that the carbon attached to the hydroxyl is not part of the CoA molecule
        match_coa = mol.GetSubstructMatch(coa_pattern)
        if match_coa is not None and hydroxy_carbon_idx in match_coa:
             continue
         
        #Validate that the hydroxyl group is in position 3 of the chain
        neighbors_alpha_carbon = [neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(alpha_carbon_idx).GetNeighbors()]
        if not thioester_carbonyl_idx in neighbors_alpha_carbon:
            continue # the chain does not connect correctly.

        neighbors_hydroxy_carbon = [neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(hydroxy_carbon_idx).GetNeighbors()]
        
        if  len([neighbor for neighbor in mol.GetAtomWithIdx(hydroxy_carbon_idx).GetNeighbors() if neighbor.GetAtomicNum() == 6]) < 2:
             continue #the hydroxy carbon should be connected to at least 2 other carbons
        if not alpha_carbon_idx in neighbors_hydroxy_carbon:
           continue # check that it is a beta hydroxy, not alpha or gamma

        # Check for a fatty acid chain (at least 3 C) attached to the carbonyl carbon of the thioester
        
        
        fatty_acid_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
        
        temp_mol = Chem.Mol(mol)
        # Remove the CoA to check the fatty acid chain
        match_coa = temp_mol.GetSubstructMatch(coa_pattern)
        if match_coa is not None:
             
             
             temp_mol = Chem.RWMol(temp_mol)
             for atom_idx in reversed(match_coa):
                  temp_mol.RemoveAtom(atom_idx)
             temp_mol = Chem.Mol(temp_mol)

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        # Remove the thioester and the 3-hydroxy group, as well as the alpha carbon of the fatty acid
        
        temp_mol_2 = Chem.RWMol(mol)
        
        
        for atom_idx in reversed([thioester_carbonyl_idx,alpha_carbon_idx, hydroxy_carbon_idx ]):
             try:
                temp_mol_2.RemoveAtom(atom_idx)
             except:
                 pass

        temp_mol_2 = Chem.Mol(temp_mol_2)
        
        
        
        fatty_matches = temp_mol_2.GetSubstructMatches(fatty_acid_chain_pattern)


        if len(fatty_matches) > 0:
             found_chain = True
             break # Found chain, can skip other matches
    
    if not found_chain:
        return False, "Missing fatty acid chain or chain too short"

    
    
    # Final validation
    return True, "Contains CoA, thioester, and a 3-hydroxy fatty acid chain."