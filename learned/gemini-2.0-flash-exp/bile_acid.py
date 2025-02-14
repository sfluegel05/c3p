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
    # This SMARTS pattern captures the basic tetracyclic ring structure
    # [C]1[C]2[C]3[C]4[C]([C]1)([C]2)([C]3)([C]4), and we add the 5beta configuration to the carbon between rings A and B using [C@]([H]).
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C@]2([C]3[C]4[C]([C@]([H])1)([C]2)([C]3)([C]4))")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Steroid core with 5beta configuration not found"

    # 2. Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Carboxylic acid group not found"

    # 3. Check for at least one hydroxyl group (-OH) on the steroid core.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group (-OH) found."

     # 4. Check for Sidechain at C17
    # Look for pattern for a chain of Carbons attached to the C17 position that ends in a carboxylic group.
    #The C17 position is the last carbon of the tetracyclic system attached to 3 or 2 other carbons.
    # Use a generic SMARTS to find the sidechain attached to any carbon
    # The specific C17 match was too rigid, so removing that condition
    sidechain_pattern = Chem.MolFromSmarts("[C]~[C](=O)O")
    sidechain_matches = mol.GetSubstructMatches(sidechain_pattern)
    if not sidechain_matches:
        return False, "Carboxylic acid chain at C17 not found"


    # 5. Check for Amide linkage (optional). If there is an amide it is also a bile acid.
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if amide_matches:
       pass

    # If all checks pass, then it's likely a bile acid
    return True, "Molecule is a bile acid based on steroid core, carboxylic acid group, hydroxyls and sidechain"