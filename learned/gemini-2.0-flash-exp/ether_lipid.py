"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid has at least one fatty acid chain attached to the glycerol backbone via an ether linkage (C-O-C).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Get the indices of the glycerol carbons
    glycerol_carbons = []
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if glycerol_matches:
         glycerol_carbons = glycerol_matches[0]
    else:
        return False, "Glycerol substructure not found"

    # Look for at least one ether group connected to a glycerol carbon (C-O-C) or (C=C-O-C), 
    # the glycerol carbon must be connected to the ether oxygen
    ether_pattern = Chem.MolFromSmarts("[CH2X4,CHX4;!$(C=O)]-[OX2]-[CX4]") # Changed to prevent ester linkages
    ether_pattern2 = Chem.MolFromSmarts("[C,c]=[C,c]-[OX2]-[CX4]") #Handles P-lipids
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    ether_matches2 = mol.GetSubstructMatches(ether_pattern2)


    has_valid_ether = False
    if ether_matches:
      for match in ether_matches:
        #check if the oxygen is connected to a glycerol carbon
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O':
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                glycerol_neighbors = [n for n in neighbors if n in glycerol_carbons]
                if glycerol_neighbors:
                    other_neighbors = [mol.GetAtomWithIdx(x) for x in neighbors if x not in glycerol_neighbors]
                    if len(other_neighbors) > 0 and other_neighbors[0].GetSymbol() == 'C':
                       has_valid_ether = True
                       break
        if has_valid_ether:
            break
    if not has_valid_ether and ether_matches2:
      for match in ether_matches2:
          #check if the oxygen is connected to a glycerol carbon
          for atom_idx in match:
             atom = mol.GetAtomWithIdx(atom_idx)
             if atom.GetSymbol() == 'O':
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                glycerol_neighbors = [n for n in neighbors if n in glycerol_carbons]
                if glycerol_neighbors:
                    other_neighbors = [mol.GetAtomWithIdx(x) for x in neighbors if x not in glycerol_neighbors]
                    if len(other_neighbors) > 0 and other_neighbors[0].GetSymbol() == 'C':
                       has_valid_ether = True
                       break
          if has_valid_ether:
            break

    if not has_valid_ether:
        return False, "No valid ether linkage found on glycerol"

    #Check for long carbon chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, f"Missing fatty acid chains"

    # Check molecular weight - lipids typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for ether lipid"

    
    return True, "Contains glycerol backbone with at least one fatty acid chain attached via ether bond"