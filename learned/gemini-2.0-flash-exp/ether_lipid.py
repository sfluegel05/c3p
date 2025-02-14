"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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
        
    # Look for at least one ether group connected to the glycerol backbone
    ether_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-O-[CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    # Look for ester groups connected to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-O-[CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    #Check if there is at least one ether linkage, but it cannot be the glycerol oxygen linked by a phosphate
    has_valid_ether = False
    for match in ether_matches:
        for atom_idx in match:
             atom = mol.GetAtomWithIdx(atom_idx)
             if atom.GetSymbol() == 'O' and atom.GetDegree() == 2:
                #check if neighbours are glycerol carbons
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                glycerol_carbons = [g.GetIdx() for g in mol.GetSubstructMatches(glycerol_pattern)[0]]
                if len(set(neighbors) & set(glycerol_carbons)) > 0 : #is it connected to glycerol?
                   #check if the next neighbour is P
                   for neighbor_idx in neighbors:
                      neighbor = mol.GetAtomWithIdx(neighbor_idx)
                      if neighbor.GetSymbol() == 'C':
                           next_neighbors = [n.GetIdx() for n in neighbor.GetNeighbors()]
                           is_phospho = False
                           for nn_idx in next_neighbors:
                              nn = mol.GetAtomWithIdx(nn_idx)
                              if nn.GetSymbol() == 'P':
                                  is_phospho = True
                           if not is_phospho:
                                has_valid_ether = True
                                break
        if has_valid_ether:
           break

    if not has_valid_ether:
      return False, "No valid ether linkage found"


    #Check for long carbon chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, f"Missing fatty acid chains"
    
    return True, "Contains glycerol backbone with at least one fatty acid chain attached via ether bond"