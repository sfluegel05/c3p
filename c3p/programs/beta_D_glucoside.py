"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
"""
Classifies: CHEBI:59826 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.
    A beta-D-glucoside is a molecule where a beta-D-glucose moiety is attached to another molecule via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-D-glucose pattern with correct stereochemistry
    beta_D_glucose_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O")
    
    # Check if the molecule contains the beta-D-glucose pattern
    if not mol.HasSubstructMatch(beta_D_glucose_pattern):
        return False, "No beta-D-glucose moiety found"

    # Get the matched atoms for the beta-D-glucose pattern
    matches = mol.GetSubstructMatches(beta_D_glucose_pattern)
    for match in matches:
        anomeric_carbon = match[0]  # C1 of the glucose moiety
        anomeric_atom = mol.GetAtomWithIdx(anomeric_carbon)
        
        # Check the glycosidic bond configuration
        for neighbor in anomeric_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen (glycosidic bond)
                # Check if the glycosidic bond is in the beta-configuration
                # In beta-D-glucosides, the glycosidic bond is axial (down)
                # We can use RDKit's stereochemistry tools to check this
                bond = mol.GetBondBetweenAtoms(anomeric_carbon, neighbor.GetIdx())
                if bond.GetStereo() == Chem.BondStereo.STEREOE:
                    # Further check the relative stereochemistry of C1 and C2
                    c2_atom = mol.GetAtomWithIdx(match[1])  # C2 of the glucose moiety
                    if anomeric_atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW and c2_atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                        return True, "Contains beta-D-glucose moiety with beta-configuration glycosidic bond"
    
    return False, "No beta-configuration glycosidic bond found"