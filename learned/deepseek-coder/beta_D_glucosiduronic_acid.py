"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid has a beta-D-glucuronic acid moiety connected via a beta-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more specific glucuronic acid pattern with stereochemistry
    glucuronic_acid_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)C(=O)O)O)O)O)[OX2]")
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid moiety found"

    # Find all matches of the glucuronic acid pattern
    matches = mol.GetSubstructMatches(glucuronic_acid_pattern)
    
    for match in matches:
        # Get the anomeric carbon (C1) and glycosidic oxygen
        anomeric_carbon = match[0]
        glycosidic_oxygen = match[5]
        
        # Get the bond between C1 and glycosidic oxygen
        bond = mol.GetBondBetweenAtoms(anomeric_carbon, glycosidic_oxygen)
        if bond is None:
            continue
            
        # Verify the beta configuration (C1-O bond is equatorial)
        # In RDKit, CHI_TETRAHEDRAL_CCW corresponds to beta configuration
        if (mol.GetAtomWithIdx(anomeric_carbon).GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW and
            bond.GetBondType() == Chem.BondType.SINGLE):
            
            # Verify the glycosidic oxygen is connected to another carbon
            glycosidic_oxygen_atom = mol.GetAtomWithIdx(glycosidic_oxygen)
            neighbors = glycosidic_oxygen_atom.GetNeighbors()
            if any(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):
                # Ensure the glycosidic oxygen is connected to another molecule
                for neighbor in neighbors:
                    if neighbor.GetIdx() != anomeric_carbon:
                        return True, "Contains beta-D-glucuronic acid moiety connected via beta-glycosidic bond"

    return False, "No beta-glycosidic bond found between glucuronic acid and another molecule"