"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: CHEBI:51868 gas molecular entity
Definition: Any main group molecular entity that is gaseous at standard temperature and
pressure (STP; 0degreeC and 100 kPa).
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get molecular weight
    mol_wt = Descriptors.MolWt(mol)
    
    # Check for elements known to be gases at STP
    gas_elements = [1, 2, 5, 6, 7, 8, 9, 10, 17, 18, 36, 54, 86]  # H, He, B, C, N, O, F, Ne, Cl, Ar, Kr, Xe, Rn
    element_nums = set([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    if len(element_nums.intersection(set(gas_elements))) == len(element_nums):
        if mol_wt < 100:
            return True, "Molecule contains only elements known to be gases at STP, and has low molecular weight"
        else:
            return False, "Molecule contains only gas elements but has high molecular weight, likely not a gas"
    
    # Check for small organic compounds
    if mol_wt < 100 and all(atom.GetAtomicNum() in [1, 5, 6, 7, 8, 9, 17] for atom in mol.GetAtoms()):
        return True, "Small organic compound, likely gaseous at STP"
    
    # Check for inorganic compounds known to be gases
    inorganic_gases = ["[O-][O+]=O", "O=C=O", "[C-]#[O+]", "ClCl", "FF", "Cl[H]", "I[H]"]
    if Chem.MolToSmiles(mol) in inorganic_gases:
        return True, "Inorganic compound known to be gaseous at STP"
    
    # If none of the above conditions are met, assume not a gas
    return False, "Not recognized as a gas molecular entity"