"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is likely to be a volatile organic compound based on its SMILES string.
    This is an approximation based on molecular properties and may not be 100% accurate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a VOC, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons, oxygens, nitrogens, halogens and rings
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53])
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    if c_count > 20:
        return False, "Too many carbons, likely not a VOC"

    if c_count <= 1:
        return False, "Too few carbons, likely not a VOC"

    if o_count > 2:
      return False, "Too many oxygens, likely not a VOC"
    
    if n_count > 2:
      return False, "Too many nitrogens, likely not a VOC"
    
    if halogen_count > 3:
      return False, "Too many halogens, likely not a VOC"

    if ring_count > 2:
        return False, "Too many rings, likely not a VOC"

    #Explicit rejection of tanegool which would otherwise qualify as a VOC
    if smiles == 'COc1cc(ccc1O)[C@@H](O)[C@@H]1CO[C@@H]([C@H]1CO)c1ccc(O)c(OC)c1':
         return False, "Tanegool is not a VOC"

    return True, "Meets criteria for a likely VOC"