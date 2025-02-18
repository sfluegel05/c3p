"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: CHEBI:137762 beta-D-galactoside
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside has a beta-configured D-galactose moiety linked via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for beta-D-galactopyranose with glycosidic bond
    # The pattern matches the galactose ring with correct stereochemistry and an O-linked group
    beta_D_galactose_pattern = Chem.MolFromSmarts("[C@@H]1(O[C@H](CO)[C@H](O)[C@@H](O)[C@H]1O)")
    if beta_D_galactose_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # Check if the molecule contains the beta-D-galactose moiety
    if not mol.HasSubstructMatch(beta_D_galactose_pattern):
        return False, "No beta-D-galactopyranose moiety found"
    
    # Additionally, verify that the anomeric oxygen is part of a glycosidic bond (connected to another atom)
    # The anomeric carbon is the one attached to the glycosidic oxygen (O in the pattern)
    # Iterate over the matches to check the glycosidic bond
    matches = mol.GetSubstructMatches(beta_D_galactose_pattern)
    for match in matches:
        anomeric_carbon_idx = match[0]  # First atom in the SMARTS is the anomeric carbon
        anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_idx)
        # Check if the anomeric carbon is bonded to an oxygen that is connected to another atom (not in the ring)
        for neighbor in anomeric_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                # Check if this oxygen is connected to a non-hydrogen atom outside the galactose ring
                for o_neighbor in neighbor.GetNeighbors():
                    if o_neighbor.GetIdx() not in match:
                        return True, "Contains beta-D-galactopyranose with glycosidic bond"
    
    return False, "No glycosidic bond found on beta-D-galactopyranose moiety"