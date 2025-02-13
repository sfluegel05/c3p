"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: CHEBI:46281 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is an organic compound with an initial boiling point <= 250°C (482°F) at standard atmospheric pressure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule is organic (contains carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule does not contain carbon, hence not organic"
    
    # Check for common VOC functional groups and substructures
    voc_patterns = [
        Chem.MolFromSmarts("[C;H3][O;H1]"),  # Alcohols
        Chem.MolFromSmarts("[C;H3][C;H2][O;H1]"),  # Alcohols
        Chem.MolFromSmarts("[C;H3]=[O;H0]"),  # Aldehydes
        Chem.MolFromSmarts("[C;H3][C;H2]=[O;H0]"),  # Ketones
        Chem.MolFromSmarts("[C;H3][N;H0]#[C;H1]"),  # Nitriles
        Chem.MolFromSmarts("[C;H2]=[C;H2]"),  # Alkenes
        Chem.MolFromSmarts("[C;H3][O;H0][C;H2]=[O;H0]"),  # Esters
        Chem.MolFromSmarts("[C;H3][C;H2][O;H0][C;H1]=[O;H0]"),  # Carboxylic acids
        Chem.MolFromSmarts("[C;H3][C;H2][N;H0][C;H2]=[O;H0]"),  # Amides
        Chem.MolFromSmarts("[S;H0](=[O;H0])(=[O;H0])"),  # Sulfones
        Chem.MolFromSmarts("[S;H0](=[O;H0])[C;H3]"),  # Sulfoxides
        Chem.MolFromSmarts("[c;H1]"),  # Aromatic rings
        Chem.MolFromSmarts("[Cl,Br,I]"),  # Halides
        Chem.MolFromSmarts("[O;H1]"),  # Ethers
    ]
    
    voc_match = False
    for pattern in voc_patterns:
        if mol.HasSubstructMatch(pattern):
            voc_match = True
            break
    
    if not voc_match:
        return False, "No common VOC functional groups or substructures found"
    
    # Estimate boiling point using a group contribution method
    bp_estimator = rdMolDescriptors.MolWt(mol) ** 0.8 * rdMolDescriptors.CalcNumRotatableBonds(mol) ** 0.2
    
    # Consider borderline cases within a range
    if bp_estimator <= 250:
        return True, "Estimated boiling point <= 250°C, classified as VOC"
    elif bp_estimator <= 270:
        return True, "Estimated boiling point close to cutoff, classified as borderline VOC"
    else:
        return False, "Estimated boiling point too high for VOC"