"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA is a CoA molecule linked to a short-chain fatty acid via a thioester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for CoA fragments
    # Flexible pyrophosphate SMARTS
    pyrophosphate_smarts = "[P](=O)([O-])[O]-[P](=O)([O-])[O-]"
    pyrophosphate_pattern = Chem.MolFromSmarts(pyrophosphate_smarts)
    pyrophosphate_matches = mol.GetSubstructMatches(pyrophosphate_pattern)


    # Ribose SMARTS - modified to also find the ribose when connected to the phosphate.
    ribose_smarts = "C1[CH](O)[CH]([CH](O)[CH](O)C1)-O-[P]"
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)


    # Pantetheine fragment SMARTS
    pantetheine_smarts = "[CX4]-C(C)(C)-[CX4]-[NX3]-[CX3](=[OX1])-[NX3]-[CX2]-[CX2]-S"
    pantetheine_pattern = Chem.MolFromSmarts(pantetheine_smarts)
    pantetheine_matches = mol.GetSubstructMatches(pantetheine_pattern)


    if not pyrophosphate_matches or not ribose_matches or not pantetheine_matches:
         return False, "Missing CoA core fragment(s)"

    # Check if fragments are connected, this could be enhanced with bond tracing
    # Get the index of the phosphate connecting the ribose and pyrophosphate
    ribose_phosphate_index = -1
    for match in ribose_matches:
        for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == "P":
                ribose_phosphate_index = atom_index
                break
        if ribose_phosphate_index != -1:
            break
    
    connected_to_pyrophosphate = False
    for match in pyrophosphate_matches:
        for atom_index in match:
             atom = mol.GetAtomWithIdx(atom_index)
             if atom.GetSymbol() == "O":
               for neighbor in atom.GetNeighbors():
                  if neighbor.GetIdx() == ribose_phosphate_index:
                      connected_to_pyrophosphate = True
                      break
        if connected_to_pyrophosphate:
              break

    if not connected_to_pyrophosphate:
      return False, "Ribose not connected to the pyrophosphate"
      
    # Define the thioester bond
    thioester_smarts = "[CX3](=[OX1])[SX2]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester bond found"

    # Check if thioester is attached to the pantetheine fragment
    attached_to_pantetheine = False
    for thioester_match in thioester_matches:
      for atom_index in thioester_match:
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetSymbol() == "S":
           for neighbor in atom.GetNeighbors():
             for match in pantetheine_matches:
                 if neighbor.GetIdx() in match:
                      attached_to_pantetheine = True
                      break
           if attached_to_pantetheine:
               break
      if attached_to_pantetheine:
              break


    if not attached_to_pantetheine:
       return False, "Thioester not attached to pantetheine fragment"



    # Find the acyl chain and count carbons.
    acyl_carbon_count = 0
    for thioester_match in thioester_matches:
        carbonyl_carbon_index = -1
        for atom_index in thioester_match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == "C" and atom.GetHybridization() == Chem.HybridizationType.SP2:
                carbonyl_carbon_index = atom_index
                break

        if carbonyl_carbon_index != -1:
             carbonyl_atom = mol.GetAtomWithIdx(carbonyl_carbon_index)
             chain_atoms = set()
             queue = [carbonyl_carbon_index]
             
             while queue:
                 current_index = queue.pop(0)
                 if current_index in chain_atoms:
                     continue
                 chain_atoms.add(current_index)
                 current_atom = mol.GetAtomWithIdx(current_index)
                 if current_atom.GetSymbol() == "C" :
                    acyl_carbon_count += 1
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetSymbol() == "C" and neighbor.GetIdx() not in chain_atoms and neighbor.GetIdx() != carbonyl_carbon_index:
                            queue.append(neighbor.GetIdx())


    # Remove the carbonyl carbon itself
    acyl_carbon_count -= 1

    if acyl_carbon_count < 1 or acyl_carbon_count > 7:
        return False, f"Fatty acid chain has {acyl_carbon_count} carbons, must be between 1 and 7"
    
    return True, "Molecule is a short-chain fatty acyl-CoA"