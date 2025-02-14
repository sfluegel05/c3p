"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate is a glucuronic acid where the carboxy group is deprotonated and with beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific SMARTS pattern for beta-D-glucuronic acid.
    # Note the use of [O;X2] for the oxygen linking to the aglycone.
    # The anomeric carbon is explictly specified as @H
    glucuronic_acid_smarts = "[O;X2][C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)C([O-])=O)O)O)O)"
    pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)

    if not mol.HasSubstructMatch(pattern):
        return False, "Molecule does not contain beta-D-glucuronic acid core"
    
    # Find the match and check anomeric carbon configuration and negative charge
    matches = mol.GetSubstructMatches(pattern)
    for match in matches:
      #Get indices
      anomeric_carbon_index = match[1] 
      ring_oxygen_index = match[0]
      carboxy_carbon_index = match[5]
      
      # Check for negative charge
      carboxy_carbon = mol.GetAtomWithIdx(carboxy_carbon_index)
      if carboxy_carbon.GetFormalCharge() != -1:
          return False, "Carboxyl group not deprotonated"


      # Get the anomeric carbon, the ring oxygen, and the carbon attached to it
      anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_index)
      ring_oxygen = mol.GetAtomWithIdx(ring_oxygen_index)

      # Get neighbor carbon attached to the anomeric carbon other than the ring oxygen
      for neighbor in anomeric_carbon.GetNeighbors():
          if neighbor.GetIdx() != ring_oxygen_index:
              neighbor_carbon_index = neighbor.GetIdx()
              break

      env = Chem.GetAtomEnvironment(mol, anomeric_carbon_index, radius=1) #Get the atoms connected to the anomeric carbon
      
      # Get the atoms in the environment (ordered)
      env_atoms = [mol.GetAtomWithIdx(idx) for idx in sorted(env.GetAtoms())] #get the atoms in the environment and order them based on index
      env_indices = [atom.GetIdx() for atom in env_atoms]
      
      # Check if the environment matches the correct stereochemistry.
      if env_indices != [ring_oxygen_index, anomeric_carbon_index, neighbor_carbon_index]: #this is a safeguard to check that indeed the correct atoms were picked
        return False, "Incorrect atom environment"

      if anomeric_carbon.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return False, "Anomeric carbon is not beta"


    return True, "Molecule is a beta-D-glucosiduronate"