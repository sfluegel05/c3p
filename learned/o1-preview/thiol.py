"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH)
    or a thioether group (-S-) is attached to a carbon atom of any
    aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to ensure accurate matching
    mol = Chem.AddHs(mol)

    # Exclude oxidized sulfur patterns (sulfoxides, sulfones, sulfonic acids, sulfates)
    oxidized_sulfur_smarts = """
        [$([#16](=O)=O),    # Sulfate S(=O)(=O)
         $([#16](=O)),      # Sulfoxide or sulfone S=O
         $([#16](=O)(=O)),  # Sulfone S(=O)(=O)
         $([#16](=O)(=O)[O-]), # Sulfonate S(=O)(=O)[O-]
         $([#16](=O)(=O)O), # Sulfonic acid S(=O)(=O)O
         $([#16](=O)(=O)N)  # Sulfonamide S(=O)(=O)N
        ]
    """
    oxidized_sulfur_pattern = Chem.MolFromSmarts(oxidized_sulfur_smarts)
    if mol.HasSubstructMatch(oxidized_sulfur_pattern):
        return False, "Contains oxidized sulfur groups (e.g., sulfoxides, sulfones), not a thiol"

    # Exclude disulfides (S-S bonds)
    disulfide_smarts = "S-S"
    disulfide_pattern = Chem.MolFromSmarts(disulfide_smarts)
    if mol.HasSubstructMatch(disulfide_pattern):
        return False, "Contains disulfide bond, not a thiol"

    # Exclude sulfur atoms bonded to heteroatoms (O, N) other than carbon and hydrogen
    unwanted_sulfur_smarts = "[#16]-[!#6;!#1]"  # Sulfur bonded to any atom that is not carbon or hydrogen
    unwanted_sulfur_pattern = Chem.MolFromSmarts(unwanted_sulfur_smarts)
    if mol.HasSubstructMatch(unwanted_sulfur_pattern):
        return False, "Contains sulfur bonded to heteroatoms other than carbon/hydrogen"

    # Define SMARTS pattern for sulfur attached to carbon atoms
    sulfur_smarts = "[#16;D2]-[#6]"  # Sulfur with two single bonds, one to carbon
    sulfur_pattern = Chem.MolFromSmarts(sulfur_smarts)

    # Search for sulfur attached to carbon atoms
    if mol.HasSubstructMatch(sulfur_pattern):
        return True, "Contains sulfur atom attached to carbon atom(s)"
    else:
        return False, "No sulfur atom attached to carbon atom found"