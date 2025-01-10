"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside must have a D-glucose moiety connected to another molecule via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for D-glucose core (pyranose form)
    # Matches the core structure with correct stereochemistry, allowing for substitutions
    glucose_core = """[C@@H]1([OR0])[C@H]([OR0])[C@@H]([OR0])[C@H]([OR0])[C@H]([C@H]1[OR0])CO[R0]"""
    glucose_pattern = Chem.MolFromSmarts(glucose_core)
    
    if not glucose_pattern or not mol.HasSubstructMatch(glucose_pattern):
        return False, "No D-glucose moiety found"

    # Find all glucose matches
    glucose_matches = mol.GetSubstructMatches(glucose_pattern)
    
    # Check each glucose match for glycosidic linkage
    for match in glucose_matches:
        anomeric_carbon = match[0]  # First atom in the pattern is the anomeric carbon
        anomeric_oxygen = None
        
        # Get the oxygen attached to anomeric carbon
        for bond in mol.GetAtomWithIdx(anomeric_carbon).GetBonds():
            other_atom = bond.GetOtherAtom(mol.GetAtomWithIdx(anomeric_carbon))
            if other_atom.GetAtomicNum() == 8:  # Oxygen
                anomeric_oxygen = other_atom
                break
                
        if anomeric_oxygen is None:
            continue
            
        # Check if oxygen is connected to something other than the glucose ring
        for bond in anomeric_oxygen.GetBonds():
            other_atom = bond.GetOtherAtom(anomeric_oxygen)
            if other_atom.GetIdx() not in match:
                # This is a glycosidic bond
                # Determine alpha/beta configuration
                chiral_tag = mol.GetAtomWithIdx(anomeric_carbon).GetChiralTag()
                if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    anomeric_config = "beta"
                else:
                    anomeric_config = "alpha"
                    
                # Determine what type of molecule it's connected to
                if other_atom.GetAtomicNum() == 6:  # Carbon
                    glycoside_type = "O-glycoside to carbon"
                elif other_atom.GetAtomicNum() == 7:  # Nitrogen
                    glycoside_type = "N-glycoside"
                elif other_atom.GetAtomicNum() == 16:  # Sulfur
                    glycoside_type = "S-glycoside"
                else:
                    glycoside_type = "O-glycoside"
                
                return True, f"Contains {anomeric_config}-D-glucose in {glycoside_type} form"

    return False, "No glycosidic linkage found"