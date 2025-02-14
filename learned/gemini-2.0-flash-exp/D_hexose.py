"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a monosaccharide with six carbons and the D configuration at carbon 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count != 6:
        return False, f"Not a hexose. Found {c_count} carbon atoms"
    if o_count < 5:
         return False, f"Not a hexose. Found {o_count} oxygen atoms"

    # Identify the carbon chain with 4 chiral centers and an alcohol at one end using SMARTS
    pattern = Chem.MolFromSmarts("[C@H](O)[C](O)[C](O)[C](O)CO")
    matches = mol.GetSubstructMatches(pattern)

    if not matches:
      pattern = Chem.MolFromSmarts("[C@@H](O)[C](O)[C](O)[C](O)CO")
      matches = mol.GetSubstructMatches(pattern)
      if not matches:
          return False, "Could not find a suitable carbon chain with 4 chiral centers."


    # Check the stereo configuration of the first chiral carbon atom in the pattern
    for match in matches:
        first_chiral_carbon_index = match[0]
        first_chiral_carbon = mol.GetAtomWithIdx(first_chiral_carbon_index)
        stereo = first_chiral_carbon.GetStereo()
        if stereo == Chem.rdchem.StereoType.STEREO_D:
          return True, "Confirmed D-configuration at carbon 5"
        
    return False, f"No carbon with the D-configuration at carbon 5 was found"