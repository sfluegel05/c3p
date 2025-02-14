"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Monosaccharide ring patterns (generic 5 or 6 membered ring with O and C)
    # This allows for substitutions on the ring
    monosaccharide_pattern1 = Chem.MolFromSmarts("[CX4,CX3][CX4,CX3][CX4,CX3][OX2][CX4,CX3][CX4,CX3]") # 6 membered ring
    monosaccharide_pattern2 = Chem.MolFromSmarts("[CX4,CX3][CX4,CX3][OX2][CX4,CX3][CX4,CX3]") # 5 membered ring
    
     # Open-chain polyhydroxy aldehyde/ketone monosaccharide pattern
    open_chain_pattern = Chem.MolFromSmarts("[CX4H](O)[CX4H](O)[CX4H](O)[CX3](=O)[CX4]")
    open_chain_pattern2 = Chem.MolFromSmarts("[CX4H](O)[CX4H](O)[CX4H](O)[CX4][CX3](=O)")
    open_chain_pattern3 = Chem.MolFromSmarts("[CX4H](O)[CX4H](O)[CX4H](O)[CX4H](O)[CX3](=O)[CX4]")
    open_chain_pattern4 = Chem.MolFromSmarts("[CX4H](O)[CX4H](O)[CX4H](O)[CX4H](O)[CX4][CX3](=O)")
    open_chain_pattern5 = Chem.MolFromSmarts("[CX4H](O)[CX4H](O)[CX4H](O)[CX4H](O)[CX4H](O)[CX3](=O)[CX4]")
    open_chain_pattern6 = Chem.MolFromSmarts("[CX4H](O)[CX4H](O)[CX4H](O)[CX4H](O)[CX4H](O)[CX4][CX3](=O)")

    
    has_monosaccharide = mol.HasSubstructMatch(monosaccharide_pattern1) or \
                        mol.HasSubstructMatch(monosaccharide_pattern2) or \
                        mol.HasSubstructMatch(open_chain_pattern) or \
                        mol.HasSubstructMatch(open_chain_pattern2) or\
                        mol.HasSubstructMatch(open_chain_pattern3) or \
                        mol.HasSubstructMatch(open_chain_pattern4) or \
                        mol.HasSubstructMatch(open_chain_pattern5) or \
                        mol.HasSubstructMatch(open_chain_pattern6)

    if not has_monosaccharide:
        return False, "No monosaccharide structure found"
    
    # Phosphate group pattern (ester) directly linked to monosaccharide
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2])([OX2])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Check for direct linkage of phosphate to a monosaccharide oxygen (ester bond)
    linked = False

    for match in phosphate_matches:
        for p_idx in match: #iterate through atoms in the phosphate substructure
            p_atom = mol.GetAtomWithIdx(p_idx)
            if p_atom.GetAtomicNum() == 8: # get only the oxygens from the phosphate
                for neighbor in p_atom.GetNeighbors():
                  if neighbor.GetAtomicNum() == 15: #check if linked to P
                    for n2 in neighbor.GetNeighbors():
                        if n2.GetAtomicNum() == 8: #check for O in P group
                            for n3 in n2.GetNeighbors():
                                if n3.GetAtomicNum() == 6:
                                    for n4 in n3.GetNeighbors():
                                      if n4.GetAtomicNum() == 8 and n4.GetDegree() == 2 and not n4.IsInRing(): #Check if it is a monosaccharide hydroxyl oxygen
                                        linked= True
                                        break
                                    if linked:
                                      break
                            if linked:
                              break
                    if linked:
                        break
                if linked:
                  break
        if linked:
            break
        
    
    if not linked:
        return False, "Phosphate not directly linked to a monosaccharide hydroxyl group via ester bond"

    return True, "Contains a monosaccharide with a phosphate group attached via ester bond"