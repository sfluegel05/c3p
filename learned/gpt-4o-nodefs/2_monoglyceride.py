"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride consists of a glycerol backbone, with a single fatty acid attached
    at the 2 position (middle) via an ester bond, and hydroxyl groups on the 1 and 3 positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a pattern for the glycerol backbone with a middle ester linkage
    glycerol_2_mono_ester_pattern = Chem.MolFromSmarts("OCC(CO)OC(=O)C")

    if mol.HasSubstructMatch(glycerol_2_mono_ester_pattern):
        # Get matches to verify the structure in detail
        matches = mol.GetSubstructMatches(glycerol_2_mono_ester_pattern)
        
        # Validate the ester group and esterification at the glycerol 2 position
        for match in matches:
            ester_carbon_idx = match[2]  # Index of the central esterified carbon
            ester_carbon = mol.GetAtomWithIdx(ester_carbon_idx)
            neighbor_atoms = [neighbor.GetSymbol() for neighbor in ester_carbon.GetNeighbors()]

            if 'O' in neighbor_atoms and len(neighbor_atoms) == 2:
                # Confirm the attachment pattern fits expected for a 2-monoglyceride
                return True, "Structure matches 2-monoglyceride with ester linkage at glycerol 2 position"
    
    return False, "Structure does not match 2-monoglyceride pattern"

# Example usage
smiles_example = "O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C(CO)CO"  # Example SMILES from dataset
result, reason = is_2_monoglyceride(smiles_example)
print(f"Classification: {result}, Reason: {reason}")