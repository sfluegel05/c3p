"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: CHEBI:33373 chalcogen
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_chalcogen(smiles: str):
    """
    Determines if a molecule is a chalcogen based on its SMILES string.
    A chalcogen is any p-block element belonging to the group 16 family of the periodic table.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcogen, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude known non-chalcogen species
    if smiles in ["[H]OO[H]", "[O-][O+]=O"]:
        return False, "The molecule is a known non-chalcogen species"
    
    # Check for chalcogen elements
    chalcogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in [8, 16, 34, 52, 84, 118]]
    if not chalcogen_atoms:
        return False, "The molecule does not contain any chalcogen elements"
    
    # Check for common chalcogen functional groups
    chalcogen_patterns = ['[OX2H]',  # Alcohols and phenols
                          '[OX2H0]',  # Ethers
                          '[SX2H]',  # Thiols
                          '[SX2H0]',  # Thioethers
                          '[SeX2H]',  # Selenols
                          '[SeX2H0]', # Selenoethers
                          '[OX1H0]',  # Oxides
                          '[SX1H0]',  # Sulfides
                          '[SeX1H0]', # Selenides
                          '[OX1H1]',  # Hydroxy groups
                          '[SX1H1]',  # Mercaptans
                          '[SeX1H1]', # Selenols
                          '[OX1H2]',  # Water
                          '[OX1H0-]', # Oxoanions
                          '[SX1H0-]', # Sulfanions
                          '[SeX1H0-]' # Selenanions
                         ]
    
    chalcogen_groups = []
    for pattern in chalcogen_patterns:
        chalcogen_groups.extend(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
    
    if chalcogen_groups:
        return True, "The molecule contains chalcogen functional groups"
    else:
        return False, "The molecule does not contain any recognized chalcogen functional groups"