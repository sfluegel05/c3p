from rdkit import Chem
from rdkit.Chem import AllChem

def is_vitamin_B1(smiles: str):
    """
    Determines if a molecule is vitamin B1 (thiamine) or a derivative.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is vitamin B1 or derivative, False otherwise
        str: Reason for classification
    """
    # Split into components in case of salts/hydrates
    components = smiles.split('.')
    
    # Find the main thiamine component
    thiamine_found = False
    thiamine_mol = None
    for component in components:
        mol = Chem.MolFromSmiles(component)
        if mol is None:
            continue

        # Look for thiazolium ring connected to pyrimidine
        thiamine_pattern = Chem.MolFromSmarts('[#6]1[#6]([#6])[S][#6][N+]1[#6][#6]2[#6]=N[#6]([#6])=N[#6]([N,NH2,NH3+])=[#6]2')
        if mol.HasSubstructMatch(thiamine_pattern):
            thiamine_found = True
            thiamine_mol = mol
            break

    if not thiamine_found:
        return False, "Not a thiamine structure"

    # Look for modifications
    modifications = []
    
    # Check for phosphate groups
    if any('P' in comp for comp in components):
        if 'OP(O)(=O)OP(O)(=O)OP' in smiles or 'OP([O-])(=O)OP([O-])(=O)OP' in smiles:
            modifications.append('triphosphate')
        elif 'OP(O)(=O)OP' in smiles or 'OP([O-])(=O)OP' in smiles:
            modifications.append('diphosphate')
        elif 'OP(O)' in smiles or 'OP([O-])' in smiles:
            modifications.append('monophosphate')

    # Check for carboxylic acid/carboxylate
    if 'C(O)=O' in smiles or 'C([O-])=O' in smiles:
        modifications.append('carboxylic acid/carboxylate')

    # Check for aldehyde 
    if '[H]C(=O)' in smiles:
        modifications.append('aldehyde')
        
    # Check for alcohol
    alcohol_pattern = Chem.MolFromSmarts('CC[OH]')
    if thiamine_mol.HasSubstructMatch(alcohol_pattern):
        modifications.append('alcohol')

    # Check for counter-ions and hydrates
    if '[Cl-]' in smiles:
        modifications.append('chloride salt')
    if '.O' in smiles:
        modifications.append('hydrate')

    if modifications:
        return True, f"Thiamine derivative with {', '.join(modifications)}"
    return True, "Thiamine (vitamin B1)"
# Pr=None
# Recall=0.0