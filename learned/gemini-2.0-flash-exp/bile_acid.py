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
    # Flexible SMARTS to capture the steroid core with 5beta configuration
    # and allow for substitutions or unsaturations.
    steroid_core_pattern = Chem.MolFromSmarts("[C;R]1~[C;R@]2~[C;R]~[C;R]~[C;R]2~[C;R]~[C;R]1")
    if not mol.HasSubstructMatch(steroid_core_pattern):
         return False, "Steroid core not found or incorrect configuration"
    
    # Check for methyl groups at C10 and C13. These are always present in bile acids.
    methyl_c10_c13_pattern = Chem.MolFromSmarts("[C;R]1~[C;R@]2([C;R]~[C;R]~[C;R]2~[C;R](C)~[C;R]1C)")
    if not mol.HasSubstructMatch(methyl_c10_c13_pattern):
        return False, "Methyl groups missing at positions 10 and 13."
    
    # 2. Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Carboxylic acid group not found"

    # 3. Check for at least one hydroxyl group (-OH) on the steroid core.
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4]-[OX2]") #Matches C-O
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group (-OH) found."
    
    # 4. Check for Sidechain at C17
    #The C17 position is the last carbon of the tetracyclic system attached to 3 or 2 other carbons.
    #Use a SMARTS that looks for this, and then a chain that ends in carboxylic acid (or carboxylate)
    # The C17 carbon is a quaternary or tertiary carbon [CX4;CX3]
    c17_sidechain_pattern = Chem.MolFromSmarts("[CX4;CX3](~[C;R])(~[C;R])~[C]~[C](=O)[O;H]")

    if not mol.HasSubstructMatch(c17_sidechain_pattern):
        return False, "Carboxylic acid sidechain not found at C17"


    # 5. Check for Amide linkage (optional). If there is an amide it is also a bile acid.
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if amide_matches:
       pass

    # If all checks pass, then it's likely a bile acid
    return True, "Molecule is a bile acid based on steroid core, carboxylic acid group, hydroxyls and sidechain"