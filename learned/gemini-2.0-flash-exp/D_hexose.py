"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


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

    # Find a chain of 6 carbons. Get indices of the longest chain.
    longest_chain_indices = []
    for path in rdMolDescriptors.GetLongestChainPaths(mol):
      if len(path) == 6:
          longest_chain_indices = path
          break
    if len(longest_chain_indices) != 6:
        return False, "Could not find a 6 carbon chain"

    # Get the 5th carbon index (index 4 in the list) from the carbon chain.
    fifth_carbon_index = longest_chain_indices[4]
    
    fifth_carbon = mol.GetAtomWithIdx(fifth_carbon_index)
    
    # Ensure the 5th carbon is indeed a carbon and is chiral.
    if fifth_carbon.GetAtomicNum() != 6:
        return False, "Fifth carbon is not carbon"
    if fifth_carbon.GetDegree() != 4:
        return False, "Fifth carbon is not tetrahedral"


    # Check the stereo configuration
    stereo = fifth_carbon.GetStereo()

    if stereo == Chem.rdchem.StereoType.STEREO_D:
          return True, "Confirmed D-configuration at carbon 5"
    else:
          return False, f"No D-configuration at carbon 5. Found {stereo}"