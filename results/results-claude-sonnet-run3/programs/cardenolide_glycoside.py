from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardenolide_glycoside(smiles: str):
    """
    Determines if a molecule is a cardenolide glycoside (cardenolide with glycosyl residues at position 3).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cardenolide glycoside, False otherwise
        str: Reason for classification
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Check for steroid core with specific stereochemistry
        steroid_core = Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~2~1~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~[#6]~[#6]~4')
        if not mol.HasSubstructMatch(steroid_core):
            return False, "No steroid core found"

        # Check for butenolide ring (unsaturated 5-membered lactone at C17)
        butenolide = Chem.MolFromSmarts('C=CC(=O)OC')
        if not mol.HasSubstructMatch(butenolide):
            return False, "No butenolide ring found"

        # Check for sugar moiety with multiple OH groups
        # This pattern looks for a pyranose ring with multiple oxygen substituents
        sugar_pattern = Chem.MolFromSmarts('O[C]1[C][C]([OH1,OR])[C]([OH1,OR])[C]([OH1,OR])[C]1')
        if not mol.HasSubstructMatch(sugar_pattern):
            return False, "No glycoside moiety found"

        # Check for glycosidic linkage at position 3
        glycosidic_linkage = Chem.MolFromSmarts('C1CCCC2(C)C1CCC3C2CC(OC4OC[C@@H]([OH1,OR])[C@H]([OH1,OR])[C@H]4[OH1,OR])CC3')
        if not mol.HasSubstructMatch(glycosidic_linkage):
            return False, "No glycosidic linkage at position 3"

        # Check for characteristic oxygen-containing substituents
        # Including 14-beta-hydroxy and other typical oxygenation patterns
        hydroxy_pattern = Chem.MolFromSmarts('[CH2]1[CH2][CH2][C@]2([OH1,OR])[C@@H]1CC[C@H]2C')
        if not mol.HasSubstructMatch(hydroxy_pattern):
            return False, "Missing characteristic hydroxyl groups"

        # Additional check for cardenolide-specific features
        cardenolide_pattern = Chem.MolFromSmarts('C1CC(=O)OC1')
        if not mol.HasSubstructMatch(cardenolide_pattern):
            return False, "Missing cardenolide-specific structural features"

        return True, "Molecule is a cardenolide glycoside with required structural features"

    except Exception as e:
        return None, f"Error analyzing molecule: {str(e)}"
# Pr=None
# Recall=0.0