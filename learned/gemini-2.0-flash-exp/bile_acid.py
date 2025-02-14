"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are characterized by a steroid core with specific hydroxyl groups and a carboxylic acid at the end of a side chain, often conjugated with glycine or taurine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the steroid core (tetracyclic ring system) with 5beta configuration.
    # Flexible SMARTS to capture the steroid core, but doesn't fix unsaturations.
    # This pattern captures the four fused rings
    steroid_core_pattern = Chem.MolFromSmarts("[C;R]1[C;R]2[C;R]3[C;R]4[C;R]1[C;R]2[C;R]34")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Steroid core not found"
    
    # 2. Check for 5beta configuration at the C5 position
    # C5 is the carbon in the first ring, shared with the second ring
    c5_beta_pattern = Chem.MolFromSmarts("[C;R]1[C@H;R]2[C;R]~[C;R]~[C;R]2[C;R]~[C;R]1")
    if not mol.HasSubstructMatch(c5_beta_pattern):
         return False, "Incorrect stereochemistry at position 5 (not 5-beta)"

    # Check for methyl groups at C10 and C13. These are always present in bile acids.
    methyl_c10_c13_pattern = Chem.MolFromSmarts("[C;R]1~[C;R@]2([C;R]~[C;R]~[C;R]2~[C;R](C)~[C;R]1C)")
    if not mol.HasSubstructMatch(methyl_c10_c13_pattern):
        return False, "Methyl groups missing at positions 10 and 13."
    
    # 3. Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Carboxylic acid group not found"

    # 4. Check for Sidechain at C17 with a terminal carboxylic acid.
    # C17 is a quaternary or tertiary carbon [CX4;CX3]
    # The sidechain can have one or more carbons before the carboxylic acid, and the carbonyl can be directly attached or separated by a carbon.
    c17_sidechain_pattern = Chem.MolFromSmarts("[CX4;CX3](~[C;R])(~[C;R])~[C]~[C](=O)[O;H]")
    if not mol.HasSubstructMatch(c17_sidechain_pattern):
        return False, "Carboxylic acid sidechain not found at C17"

    # Check for hydroxyl groups at positions 3, 7 and/or 12
    # Note: It's not mandatory that all three are present, so they are optional
    # The -OH must be attached to a stereocenter.
    hydroxyl_3_pattern = Chem.MolFromSmarts("[C@H;R](O)~[C;R]~[C;R]")
    hydroxyl_7_pattern = Chem.MolFromSmarts("[C@H;R](O)~[C;R]~[C;R]~[C;R]~[C;R]")
    hydroxyl_12_pattern = Chem.MolFromSmarts("[C@H;R](O)~[C;R]~[C;R]~[C;R]~[C;R]~[C;R]")
    
    if mol.HasSubstructMatch(hydroxyl_3_pattern) or mol.HasSubstructMatch(hydroxyl_7_pattern) or mol.HasSubstructMatch(hydroxyl_12_pattern):
        pass #At least one hydroxyl in the typical positions is present.
    else:
        pass #It is still a bile acid if it does not have an hydroxyl at those positions.
    
    # 5. Amide linkages are optional. Do not check.

    # If all checks pass, then it's likely a bile acid
    return True, "Molecule is a bile acid based on steroid core, 5-beta configuration, carboxylic acid group, hydroxyls and sidechain"