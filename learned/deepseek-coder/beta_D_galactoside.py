"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: CHEBI:18133 beta-D-galactoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside is a galactose ring in the beta-configuration, where the hydroxyl group at the anomeric carbon (C1) is in the axial position,
    and the anomeric carbon is part of a glycosidic bond (connected to another group via an oxygen atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible beta-D-galactoside pattern
    # The pattern matches a galactose ring with the beta-configuration at the anomeric center (C1)
    # and ensures that the anomeric carbon is part of a glycosidic bond (connected to another group via an oxygen atom).
    beta_D_galactoside_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)[OX2]")

    # Check if the molecule contains the beta-D-galactoside pattern
    if mol.HasSubstructMatch(beta_D_galactoside_pattern):
        # Further verify that the anomeric carbon is connected to another group via an oxygen atom
        matches = mol.GetSubstructMatches(beta_D_galactoside_pattern)
        for match in matches:
            anomeric_carbon = match[0]
            for neighbor in mol.GetAtomWithIdx(anomeric_carbon).GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                    return True, "Contains a galactose ring in the beta-configuration with a glycosidic bond at the anomeric carbon"
        return False, "No glycosidic bond found at the anomeric carbon"
    else:
        return False, "No galactose ring in the beta-configuration found"

# Example usage:
# smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"  # Methyl beta-D-galactoside
# print(is_beta_D_galactoside(smiles))