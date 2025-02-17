"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: Hydroxy Fatty Acid
Definition: Any fatty acid carrying one or more hydroxy substituents.
A molecule qualifies if it has a carboxylic acid group and at least one additional hydroxyl group.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    
    Requirements:
      - Must contain a carboxylic acid group (C(=O)[OX2H]) to be considered a fatty acid.
      - Must contain at least one additional hydroxyl substituent (–OH) outside the acid group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a hydroxy fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to accurately detect –OH groups.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for the carboxylic acid group (i.e., fatty acid carboxyl group):
    # This matches a carbonyl carbon (C(=O)) with an –OH group.
    acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found (not a fatty acid)"
    
    # Record indices for the oxygen atoms that belong to the carboxylic acid group.
    acid_oh_indices = set()
    for match in acid_matches:
        # The pattern "C(=O)[OX2H]" returns matches where match[1] is the oxygen in the acid group.
        if len(match) >= 2:
            acid_oh_indices.add(match[1])
    
    # Define a general SMARTS pattern to detect hydroxyl (–OH) groups.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    # Filter out hydroxyls that are part of the acid group.
    non_acid_oh = [match[0] for match in oh_matches if match[0] not in acid_oh_indices]
    if not non_acid_oh:
        return False, "No hydroxy substituent found outside the carboxylic acid group"
    
    # Optional: Ensure there is a minimal number of carbons to support the fatty acid definition.
    # A simple threshold is used here (e.g., at least 4 carbons).
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) < 4:
        return False, "Too few carbons to be considered a fatty acid"
    
    # If all conditions are met, classify as a hydroxy fatty acid.
    return True, "Molecule contains a carboxylic acid group and additional hydroxy substituent(s), classifying it as a hydroxy fatty acid"