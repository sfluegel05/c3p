"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is a fatty acid with the carboxyl group esterified with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the methyl ester group pattern (C(=O)OC)
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])OC")
    methyl_ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    
    if len(methyl_ester_matches) == 0:
        return False, "No methyl ester group found"
    
    if len(methyl_ester_matches) > 1:
        return False, "Multiple methyl ester groups found, only 1 is allowed"
    
    # Check for the methyl group attached to the oxygen on the ester
    for match in methyl_ester_matches:
       ester_o_idx = -1
       for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == 'O':
                ester_o_idx = idx
                break
       if ester_o_idx != -1:
        methyl_attached = False
        for neighbor in mol.GetAtomWithIdx(ester_o_idx).GetNeighbors():
             if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() == 3:
                methyl_attached = True
                break
        if not methyl_attached:
            return False, "Methyl group not found on ester oxygen"

    # Check for fatty acid chains (long carbon chains attached to carbonyl of ester)
    fatty_acid_chain_found = False
    for match in methyl_ester_matches:
         carbonyl_idx = -1
         for idx in match:
              atom = mol.GetAtomWithIdx(idx)
              if atom.GetSymbol() == 'C' and atom.GetTotalValence() == 3:
                carbonyl_idx = idx
                break
         if carbonyl_idx != -1:
             
            alkyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
            
            # get all neighbours for carbonyl, not including ester oxygen
            neighbours = [neighbour for neighbour in mol.GetAtomWithIdx(carbonyl_idx).GetNeighbors() if mol.GetAtomWithIdx(neighbour.GetIdx()).GetSymbol() != 'O' ]
            
            
            for neighbour in neighbours:
                submol = Chem.PathToSubmol(mol, [carbonyl_idx, neighbour.GetIdx()])
                
                if submol.HasSubstructMatch(alkyl_chain_pattern):
                    fatty_acid_chain_found = True
                    break
            
            if fatty_acid_chain_found:
                break
    
    if not fatty_acid_chain_found:
          return False, "Missing fatty acid chain"


    return True, "Contains a methyl ester group and a fatty acid chain"