from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole

def is_dioxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a dioxo monocarboxylic acid anion.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a dioxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for carboxylate anion group C([O-])=O
    carboxylate_pattern = Chem.MolFromSmarts('C([O-])=O')
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate anion group found"
    
    if len(mol.GetSubstructMatches(carboxylate_pattern)) > 1:
        return False, "More than one carboxylate anion group found"

    # Check for ketone groups C(=O)
    ketone_pattern = Chem.MolFromSmarts('C(=O)C')
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    if len(ketone_matches) < 2:
        return False, "Less than two ketone groups found"
    
    # Verify ketones are not part of carboxylate
    carboxylate_match = mol.GetSubstructMatch(carboxylate_pattern)
    if carboxylate_match:
        carboxylate_carbon = carboxylate_match[0]
        
        # Filter out ketone matches that include the carboxylate carbon
        real_ketones = [match for match in ketone_matches 
                       if carboxylate_carbon not in match]
        
        if len(real_ketones) < 2:
            return False, "Less than two ketone groups found (excluding carboxylate)"
            
        return True, "Contains one carboxylate anion and two ketone groups"

    return False, "Structure does not match dioxo monocarboxylic acid anion pattern"
# Pr=0.8823529411764706
# Recall=1.0