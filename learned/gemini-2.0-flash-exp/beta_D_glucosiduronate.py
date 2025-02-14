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
    
    # Define a SMARTS pattern for the beta-D-glucuronic acid core, with beta anomeric carbon and deprotonated carboxy.
    # explicitly match the ring, also specify that the oxygen is OX2, and that it links to another atom.
    glucuronic_acid_smarts = "[OX2;H0][C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)C([O-])=O)O)O)O)-[O;H0]"

    pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)

    if not mol.HasSubstructMatch(pattern):
        return False, "Molecule does not contain beta-D-glucuronic acid core"

    # check that the anomeric carbon is indeed beta
    match = mol.GetSubstructMatch(pattern)
    if match:
      anomeric_carbon_index = match[1] #index of the anomeric carbon
      anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_index)
      if anomeric_carbon.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
          return False, "Anomeric carbon not in beta configuration"
      # Check that the anomeric oxygen is connected with a single bond
      anomeric_oxygen_index = match[0]
      anomeric_oxygen = mol.GetAtomWithIdx(anomeric_oxygen_index)
      if anomeric_oxygen.GetDegree() != 2:
         return False, "Anomeric oxygen not linked via glycosidic bond"

      # check that it is not a sulfate
      for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16 and atom.GetNeighbors() and any(neighbor.GetAtomicNum() == 8 for neighbor in atom.GetNeighbors()):
          return False, "Molecule contains a sulfate group"

    return True, "Molecule is a beta-D-glucosiduronate"