from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_citrate_salt(smiles: str):
    """
    Determines if a molecule is a citrate salt.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a citrate salt, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Split into components
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(fragments) < 2:
        return False, "Not a salt - single component"
        
    # Define SMARTS patterns for citrate anion
    citrate_patterns = [
        # Fully protonated citric acid
        'OC(=O)CC(O)(CC(O)=O)C(O)=O',
        # Triply deprotonated citrate
        '[O-]C(=O)CC(O)(CC([O-])=O)C([O-])=O',
        # Mix of protonated/deprotonated forms
        'OC(=O)C(O)(CC([O-])=O)CC([O-])=O',
        # Alternative representation
        'C(=O)(O)C(CC(=O)O)(CC(=O)O)O'
    ]
    
    # Check if any fragment matches citrate pattern
    found_citrate = False
    for fragment in fragments:
        frag_smiles = Chem.MolToSmiles(fragment)
        if any(frag_smiles == pattern for pattern in citrate_patterns):
            found_citrate = True
            break
            
    if not found_citrate:
        # Try substructure matching as backup
        for fragment in fragments:
            for pattern in citrate_patterns:
                pattern_mol = Chem.MolFromSmiles(pattern)
                if pattern_mol is None:
                    continue
                if fragment.HasSubstructMatch(pattern_mol):
                    found_citrate = True
                    break
            if found_citrate:
                break
                
    if not found_citrate:
        return False, "No citrate anion found"
        
    # Look for organic cation component
    has_organic_cation = False
    for fragment in fragments:
        if not fragment.HasSubstructMatch(Chem.MolFromSmiles(citrate_patterns[0])):
            # Check for organic atoms besides H
            organic_atoms = ['C', 'N']
            if any(atom.GetSymbol() in organic_atoms for atom in fragment.GetAtoms()):
                has_organic_cation = True
                break
                
    # Look for inorganic cations
    has_inorganic_cation = False
    inorganic_cations = ['[H+]', '[NH4+]']
    for fragment in fragments:
        frag_smiles = Chem.MolToSmiles(fragment)
        if frag_smiles in inorganic_cations:
            has_inorganic_cation = True
            break
            
    if not (has_organic_cation or has_inorganic_cation):
        return False, "No cation component found"
        
    return True, "Valid citrate salt"
# Pr=1.0
# Recall=1.0