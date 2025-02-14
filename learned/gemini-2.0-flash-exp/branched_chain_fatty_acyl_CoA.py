"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA is a fatty acyl-CoA derived from a branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for CoA
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
          return False, "CoA moiety not found"

    # Define SMARTS pattern for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, "Incorrect number of thioester linkages"

    # Find carbonyl carbon of thioester
    carbonyl_carbon_indices = [match[0] for match in mol.GetSubstructMatches(thioester_pattern)]
    if not carbonyl_carbon_indices:
       return False, "Carbonyl group of thioester not found"
    
    carbonyl_carbon_index = carbonyl_carbon_indices[0]
    
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_index)
    
    # Get neighbors of carbonyl carbon
    neighbors = [atom.GetIdx() for atom in carbonyl_carbon.GetNeighbors()]
    
    # Find the carbon that is part of the fatty acyl chain
    fatty_chain_carbon_index = None
    for neighbor_index in neighbors:
      neighbor = mol.GetAtomWithIdx(neighbor_index)
      if neighbor.GetSymbol() == 'C':
        fatty_chain_carbon_index = neighbor_index
        break
    
    if fatty_chain_carbon_index is None:
      return False, "Fatty chain carbon not found"


    # BFS to traverse the fatty acyl chain and find branching point
    queue = [fatty_chain_carbon_index]
    visited_atoms = set()
    acyl_chain_atoms = set()
    is_branched = False

    while queue:
        current_index = queue.pop(0)

        if current_index in visited_atoms:
            continue
        visited_atoms.add(current_index)

        current_atom = mol.GetAtomWithIdx(current_index)
        acyl_chain_atoms.add(current_index)

        carbon_neighbors = [neighbor.GetIdx() for neighbor in current_atom.GetNeighbors() if neighbor.GetSymbol() == 'C']
        
        # Check for branching only on the chain.
        branching_neighbors = 0
        for neighbor_index in carbon_neighbors:
            if neighbor_index not in visited_atoms and neighbor_index not in acyl_chain_atoms:
                branching_neighbors+=1

        if branching_neighbors > 1:
             is_branched = True
             break

        for neighbor_index in carbon_neighbors:
            if neighbor_index not in visited_atoms and neighbor_index not in acyl_chain_atoms:
                queue.append(neighbor_index)


    if not is_branched:
          return False, "No branch found on the fatty acyl chain."

    # Chain length check
    if len(acyl_chain_atoms) < 4:
      return False, "Chain too short to be a fatty acid."

    #Check for rings.
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        for ring in ring_info.AtomRings():
           if any(atom_index in acyl_chain_atoms for atom_index in ring):
                return False, "Acyl chain is part of a ring"

    return True, "Molecule is a branched-chain fatty acyl-CoA"