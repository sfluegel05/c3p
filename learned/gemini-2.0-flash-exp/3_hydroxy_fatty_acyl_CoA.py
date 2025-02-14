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

    # Combined pattern for thioester and 3-hydroxy fatty acid chain
    # This pattern includes the thioester, the alpha carbon, the beta carbon, and the 3-hydroxy group.
    hydroxy_fatty_acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])S[CX4][CX4]([OX2])[CX3](=[OX1])")
    
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_fatty_acyl_pattern)
    if len(hydroxy_matches) == 0:
        return False, "No 3-hydroxy fatty acyl thioester group found"
    
    # Fatty acid chain pattern attached to the 3-hydroxy group. This pattern will match carbon chains, including double bonds
    fatty_acid_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    
    
    found_chain = False
    for match in hydroxy_matches:
        # Get the atom indices that match the pattern
        thioester_carbon_idx = match[0]  # C=O of thioester
        alpha_carbon_idx = match[2]      # carbon attached to S
        hydroxy_carbon_idx = match[3] # the carbon with OH
        carbonyl_carbon_idx = match[4] # the carbon of the carbonyl

        # Check that there is a long chain linked to the hydroxy carbon
        # Build a new molecule without the CoA moiety
        
        
        
        
        
        temp_mol = Chem.RWMol(mol)
        
        # Get the atoms of the CoA
        match_coa = mol.GetSubstructMatch(coa_pattern)
        
        # Remove atoms of CoA to have better counting
        if match_coa is not None:
          for atom_idx in reversed(match_coa):
           temp_mol.RemoveAtom(atom_idx)

        
        #Now remove the thioester and 3-hydroxy group
        for atom_idx in reversed([thioester_carbon_idx, temp_mol.GetAtomWithIdx(alpha_carbon_idx).GetIdx(), temp_mol.GetAtomWithIdx(hydroxy_carbon_idx).GetIdx(), temp_mol.GetAtomWithIdx(carbonyl_carbon_idx).GetIdx()]):
             try:
                temp_mol.RemoveAtom(atom_idx)
             except:
                pass
        
        
        #Search for fatty acid chains on the remaining mol
        
        if temp_mol.GetNumAtoms() > 0:
             
             fatty_matches = temp_mol.GetSubstructMatches(fatty_acid_chain_pattern)
             
             if len(fatty_matches) >= 1:
                 
                 found_chain = True
                 break # No need to check other matches
    
    if not found_chain:
         return False, "Missing or short fatty acid chain."

    #Check for carbons on the fatty acid. Needs to be done on the fatty acid without the CoA
    c_count = 0
    if found_chain:
        
      
      
        temp_mol = Chem.RWMol(mol)
        
        # Get the atoms of the CoA
        match_coa = mol.GetSubstructMatch(coa_pattern)
        
        # Remove atoms of CoA to have better counting
        if match_coa is not None:
          for atom_idx in reversed(match_coa):
           temp_mol.RemoveAtom(atom_idx)
        
        
        #Now remove the thioester and 3-hydroxy group
        for match in hydroxy_matches:
            thioester_carbon_idx = match[0]  # C=O of thioester
            alpha_carbon_idx = match[2]      # carbon attached to S
            hydroxy_carbon_idx = match[3] # the carbon with OH
            carbonyl_carbon_idx = match[4] # the carbon of the carbonyl
            
            for atom_idx in reversed([thioester_carbon_idx, temp_mol.GetAtomWithIdx(alpha_carbon_idx).GetIdx(), temp_mol.GetAtomWithIdx(hydroxy_carbon_idx).GetIdx(), temp_mol.GetAtomWithIdx(carbonyl_carbon_idx).GetIdx()]):
               try:
                temp_mol.RemoveAtom(atom_idx)
               except:
                   pass
           
        c_count = sum(1 for atom in temp_mol.GetAtoms() if atom.GetAtomicNum() == 6)


    if c_count < 3:
       return False, "Too few carbons in the fatty acid chain."

    return True, "Contains CoA, thioester, and a 3-hydroxy fatty acid chain"