"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    A tocol is a chromanol with a chroman-6-ol skeleton substituted at position 2 by a
    saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the chroman-6-ol core, explicitly capturing the core C atom.
    # This pattern uses generic carbons ([C]) and * for any atom and ensures a quaternary carbon at position 2.
    # Includes oxygen single or double bond [OX1,OX2]
    chromanol_pattern = Chem.MolFromSmarts('[c]1[c]([OX1,OX2])[c]([c])[c]2[O][C][C][c]2[c]1-[C](*)(*)-[#6]') 

    core_matches = mol.GetSubstructMatches(chromanol_pattern)
    if not core_matches:
        return False, "Chroman-6-ol core not found"
    
    # Get the 2-position of the core.
    # Position 2 is the last matched atom, as the SMARTS pattern is built to capture it explicitly
    position_2_carbon = mol.GetAtomWithIdx(core_matches[0][-1])

    # Check for the isoprenoid chain and its length
    isoprenoid_chain_carbons = []
    visited_atoms = set()
    
    def traverse_chain(current_atom, chain):
      
      if current_atom.GetIdx() in visited_atoms:
          return
      
      visited_atoms.add(current_atom.GetIdx())
      
      if current_atom.GetSymbol() == "C":
        chain.append(current_atom)
      
      neighbors = current_atom.GetNeighbors()

      for neighbor in neighbors:
         if neighbor.GetSymbol() == "C":
           
            num_carbon_neighbors = 0
            for neighbor_of_neighbor in neighbor.GetNeighbors():
                if neighbor_of_neighbor.GetSymbol() == 'C':
                    num_carbon_neighbors += 1
            
            if num_carbon_neighbors <=3:
               traverse_chain(neighbor, chain)

    traverse_chain(position_2_carbon, isoprenoid_chain_carbons)


    isoprenoid_count = len(isoprenoid_chain_carbons)

    if isoprenoid_count != 15:
        return False, f"Isoprenoid chain does not have 15 carbons. Found {isoprenoid_count} carbon atoms in the chain."

    #Check the molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300: # Tocols are generally above 300 Da
      return False, "Molecular weight too low for tocol"


    return True, "Contains chroman-6-ol core with a 15-carbon isoprenoid chain at position 2."