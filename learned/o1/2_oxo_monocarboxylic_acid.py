"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: CHEBI:37727 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid is any monocarboxylic acid having a 2-oxo substituent.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for monocarboxylic acid group (protonated or deprotonated)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]')

    # Find matches for the carboxylic acid group
    acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(acid_matches) == 0:
        return False, "No carboxylic acid groups found"

    # Initialize flag to determine if a 2-oxo substituent is present
    found_2_oxo = False

    # Iterate over each carboxylic acid group
    for match in acid_matches:
        acid_carbon_idx = match[0]
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)

        # Get alpha carbons (carbons adjacent to the carboxylic acid carbon)
        alpha_carbons = []
        for neighbor in acid_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                alpha_carbons.append(neighbor)
        
        # Check each alpha carbon for the presence of a ketone (=O) group
        for alpha_carbon in alpha_carbons:
            # Iterate over bonds connected to the alpha carbon
            has_ketone = False
            for bond in alpha_carbon.GetBonds():
                neighbor = bond.GetOtherAtom(alpha_carbon)
                # Check for a double bond to an oxygen atom (ketone group)
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
                    has_ketone = True
                    break
            if has_ketone:
                found_2_oxo = True
                break  # Stop checking other alpha carbons
        if found_2_oxo:
            break  # Stop checking other carboxylic acid groups

    if not found_2_oxo:
        return False, "No 2-oxo substituent found on alpha carbon(s)"

    return True, "Molecule is a 2-oxo monocarboxylic acid"