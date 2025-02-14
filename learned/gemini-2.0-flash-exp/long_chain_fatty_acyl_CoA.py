"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA is a coenzyme A derivative with a fatty acid chain of 13-22 carbons,
    linked via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    #  Check for CoA substructure - using a more complete SMARTS
    coa_pattern = Chem.MolFromSmarts("N[c]1[n][c]([n][c]2[n]1[C]([C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OCC(C)(C)[C@H](O)CC(=O)NCCS)=O)O)O)=N2")
    if not mol.HasSubstructMatch(coa_pattern):
       return False, "Molecule does not contain the required coenzyme A substructure"

    # Check for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("[SX2][CX3](=[OX1])")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
         return False, "Molecule does not contain exactly one thioester group"

    # Find the carbonyl carbon of the thioester group
    carbonyl_carbon = None
    for match in thioester_matches:
        for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP2:
                carbonyl_carbon = atom
                break
        if carbonyl_carbon:
            break
            
    if carbonyl_carbon is None:
        return False, "Could not identify the carbonyl carbon"
    
    # Find the carbon atom directly attached to the carbonyl carbon (excluding the sulfur)
    fatty_acid_start = None
    for neighbor in carbonyl_carbon.GetNeighbors():
         if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in [atom.GetIdx() for match in mol.GetSubstructMatches(Chem.MolFromSmarts("[SX2]")) for atom in match]:
             fatty_acid_start = neighbor
             break

    if fatty_acid_start is None:
        return False, "Fatty acid chain not directly linked to the thioester group"


    # Trace and count carbons in the fatty acid chain using DFS and exclude Coenzyme A atoms.
    visited_atoms = set()
    stack = [fatty_acid_start]
    carbon_count = 0
    
    coa_atoms = set()
    for match in mol.GetSubstructMatches(coa_pattern):
        for idx in match:
            coa_atoms.add(idx)
    
    while stack:
      current_atom = stack.pop()
      if current_atom.GetIdx() in visited_atoms or current_atom.GetIdx() in coa_atoms:
        continue
      
      visited_atoms.add(current_atom.GetIdx())
      
      if current_atom.GetAtomicNum() == 6: # Carbon
        carbon_count+=1

      for neighbor in current_atom.GetNeighbors():
        stack.append(neighbor)
    
    # Check chain length
    if carbon_count < 13 or carbon_count > 22:
      return False, f"Fatty acid chain has {carbon_count} carbons, should be between 13 and 22."

    return True, "Molecule is a long-chain fatty acyl-CoA"