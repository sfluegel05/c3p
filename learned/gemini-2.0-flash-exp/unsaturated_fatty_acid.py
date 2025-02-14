"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains a carboxylic acid group at the end of a linear carbon chain,
    at least one C=C or C#C bond, and no aromatic rings or phosphorus atoms.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an unsaturated fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for phosphorus atoms
    has_phosphorus = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    if has_phosphorus:
        return False, "Contains phosphorus, not a fatty acid"

    # Check for aromatic rings
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(aromatic_pattern):
        return False, "Contains aromatic rings, not a fatty acid"

    # More specific check for carboxylic acid group at the end of a chain, both protonated and deprotonated, and with no explicit Hs.
    acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0-,OX1-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Must have exactly one carboxylic acid group, found {len(acid_matches)}"
    
    # Check that the chain is connected to the carboxyl group and has at least 4 carbons, including the carbonyl carbon
    chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0-,OX1-]-[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[!R;CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
      return False, "Does not contain a sufficiently long linear carbon chain attached to the carboxylic acid"

    # Check for at least one C=C or C#C bond *in the chain*
    unsaturated_chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0-,OX1-]-[CX4,CX3]~[CX4,CX3]=[CX4,CX3]~[!R;CX4,CX3]")
    unsaturated_chain_pattern2 = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0-,OX1-]-[CX4,CX3]~[CX4,CX3]#[CX4,CX3]~[!R;CX4,CX3]")
    
    if not (mol.HasSubstructMatch(unsaturated_chain_pattern) or mol.HasSubstructMatch(unsaturated_chain_pattern2)):
       return False, "Does not contain any C=C or C#C bond in the main carbon chain"


    # Check number of carbons in total (at least 4)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 4:
        return False, "Too few carbon atoms"

    # If all checks pass
    return True, "Contains a carboxylic acid group, a long chain and at least one C=C or C#C bond"