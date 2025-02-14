"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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

    # Define SMARTS for branched alkyl chain, note that we don't want a linear chain
    branched_alkyl_pattern = Chem.MolFromSmarts("[CX4]([CX4])([CX4])")
    
    # Find alkyl chain attached to carbonyl of thioester
    carbonyl_carbon_indices = [match[0] for match in mol.GetSubstructMatches(thioester_pattern)]
    if not carbonyl_carbon_indices:
       return False, "Carbonyl group of thioester not found"
    
    carbonyl_carbon_index = carbonyl_carbon_indices[0]
    
    #Find atoms attached to the carbonyl carbon of the thioester
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

    
    #Check if the chain is branched
    fatty_chain_carbon = mol.GetAtomWithIdx(fatty_chain_carbon_index)
    
    has_branch = False
    
    for atom in mol.GetAtoms():
      if atom.GetIdx() == fatty_chain_carbon_index:
        continue
      if atom.GetSymbol() == 'C' and atom.GetDegree() > 2:
        neighbors_of_atom = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol()=='C']
        if len(neighbors_of_atom) > 2:
           has_branch=True
           break
    if not has_branch:
        return False, "No branch found on the fatty acyl chain."

    #Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5: #Minimum for branched
        return False, "Too few carbons for a fatty acid chain."

    #Check molecular weight > 700
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
      return False, "Molecular weight too low for branched-chain fatty acyl-CoA"
    
    return True, "Molecule is a branched-chain fatty acyl-CoA"