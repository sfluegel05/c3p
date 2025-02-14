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
    pyrophosphate_smarts = "[P](=O)-[O]-[P](=O)"
    pyrophosphate_pattern = Chem.MolFromSmarts(pyrophosphate_smarts)
    pyrophosphate_matches = mol.GetSubstructMatches(pyrophosphate_pattern)
    if not pyrophosphate_matches:
        return False, "Missing pyrophosphate fragment"


    # Ribose SMARTS
    ribose_smarts = "[C]1[CH](O)[CH]([CH](O)[CH](O)[C]1)-[O]-[P]"
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)
    if not ribose_matches:
        return False, "Missing ribose fragment"


    # Pantetheine fragment SMARTS
    pantetheine_smarts = "[CX4]-C(C)(C)-[CX4]-[NX3]-[CX3](=[OX1])-[NX3]-[CX2]-[CX2]-S"
    pantetheine_pattern = Chem.MolFromSmarts(pantetheine_smarts)
    pantetheine_matches = mol.GetSubstructMatches(pantetheine_pattern)
    if not pantetheine_matches:
       return False, "Missing pantetheine fragment"

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

    # Check if ribose and pyrophosphate are connected
    connected_to_pyrophosphate = False
    for ribose_match in ribose_matches:
       for atom_index in ribose_match:
           atom = mol.GetAtomWithIdx(atom_index)
           if atom.GetSymbol() == "P":
              for neighbor in atom.GetNeighbors():
                 for pyrophosphate_match in pyrophosphate_matches:
                     if neighbor.GetIdx() in pyrophosphate_match:
                        connected_to_pyrophosphate = True
                        break
              if connected_to_pyrophosphate:
                   break
       if connected_to_pyrophosphate:
            break

    if not connected_to_pyrophosphate:
        return False, "Ribose not connected to pyrophosphate"


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
            # count non-hydrogen neighbors of carbonyl carbon
            acyl_carbon_count = sum(1 for neighbor in carbonyl_atom.GetNeighbors() if neighbor.GetSymbol() != "H")
            # remove the carbonyl oxygen
            acyl_carbon_count -= 1

    if acyl_carbon_count < 1 or acyl_carbon_count > 7:
        return False, f"Fatty acid chain has {acyl_carbon_count} carbons, must be between 1 and 7"
    
    return True, "Molecule is a short-chain fatty acyl-CoA"