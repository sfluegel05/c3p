"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import GetLongestChainPath

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

    # Find the longest chain of carbon atoms
    try:
      chain = GetLongestChainPath(mol)
      if len(chain) < 6:
           return False, "Carbon chain too short"
      carbon_chain = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() in chain and atom.GetAtomicNum() == 6]
      if len(carbon_chain) != 6:
          return False, "Longest chain is not of length 6"
    except Exception:
      return False, "Could not identify a carbon chain."
    
    # Check for D configuration at C5
    try:
         c5_index = carbon_chain[4] # Index 4 corresponds to the 5th carbon
         c5_atom = mol.GetAtomWithIdx(c5_index)
         stereo = c5_atom.GetStereo()
         
         if stereo == Chem.rdchem.StereoType.STEREO_D:
             return True, "CORRECT D-hexose configuration confirmed"
         else:
              return False, f"C5 configuration is not D, but is {stereo}"
    except Exception as e:
         return False, f"Error checking chirality at C5: {e}"