"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    A 1-O-acylglycerophosphoethanolamine has a glycerol backbone with a phosphate group
    attached to the 3-position, esterified with ethanolamine and a fatty acid at the 1-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for the glycerol backbone with phosphate and ethanolamine at 3-position
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4]([OX2])[CHX4](O)[CH2X4]OP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with phosphate and ethanolamine at 3-position not found"
    
    # Find the attachment point on glycerol (oxygen attached to C1 of glycerol)
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "Could not identify attachment point of fatty acid"
    
    glycerol_match = matches[0]
    c1_glycerol = glycerol_match[0] #The first carbon of glycerol, it's the C connected to the ester oxygen
    
    #Find the ester oxygen connected to C1 of glycerol
    ester_oxygen = None
    for neighbor in mol.GetAtomWithIdx(c1_glycerol).GetNeighbors():
      if neighbor.GetSymbol() == 'O':
         ester_oxygen = neighbor
         break
    if ester_oxygen is None:
        return False, "Could not find ester oxygen at C1 of glycerol"

    # Find the fatty acid chain connected to the ester oxygen.
    # We use a SMARTS to find an acyl group which is C(=O)- chain
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[#6]")
    acyl_matches = []
    for neighbor in ester_oxygen.GetNeighbors():
        temp_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        #Create a new molecule and select only the neighbor.
        temp_mol_atoms = [atom for atom in temp_mol.GetAtoms()]
        to_remove = [atom for atom in temp_mol_atoms if atom.GetIdx() != neighbor.GetIdx()]
        for atom in to_remove:
           temp_mol.RemoveAtom(atom.GetIdx())
        
        if temp_mol.HasSubstructMatch(acyl_pattern):
            acyl_matches.append(neighbor.GetIdx())
    if len(acyl_matches) != 1:
        return False, "Could not find an acyl group connected to the ester oxygen"

    acyl_carbon = acyl_matches[0]

    # Check that this carbon is connected to the ester oxygen
    if not mol.GetAtomWithIdx(acyl_carbon).HasBondWith(ester_oxygen):
        return False, "Acyl carbon not connected to the ester oxygen"

    # Get the indices of atoms of the acyl chain
    acyl_chain_atoms = [acyl_carbon]
    queue = [acyl_carbon]
    visited = {acyl_carbon}

    while queue:
      current_atom_idx = queue.pop(0)
      current_atom = mol.GetAtomWithIdx(current_atom_idx)
      for neighbor in current_atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx not in visited and neighbor.GetSymbol() != 'O': #only follow C chain, not other O
          visited.add(neighbor_idx)
          acyl_chain_atoms.append(neighbor_idx)
          queue.append(neighbor_idx)
          
    # Count rotatable bonds within the fatty acid chain
    n_rotatable = 0
    for bond in mol.GetBonds():
       if bond.GetBeginAtomIdx() in acyl_chain_atoms and bond.GetEndAtomIdx() in acyl_chain_atoms and bond.IsRotatable():
        n_rotatable +=1
        
    if n_rotatable < 4:
      return False, "Fatty acid chain at the 1-position too short"

    #Check for at least one phosphorus
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
       return False, "Must have exactly one phosphorus"


    return True, "Contains glycerol backbone with phosphate and ethanolamine at the 3-position and a fatty acid chain at the 1-position"