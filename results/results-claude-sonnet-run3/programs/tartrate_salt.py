from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops

def is_tartrate_salt(smiles: str):
    """
    Determines if a molecule is a tartrate salt.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a tartrate salt, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Split into fragments
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    
    # Look for tartrate anion fragment
    tartrate_patterns = [
        # Tartrate patterns with different protonation states
        'O[C@H]([C@@H](O)C(O)=O)C(O)=O',  # Tartaric acid
        'O[C@H]([C@@H](O)C([O-])=O)C([O-])=O', # Tartrate dianion
        'O[C@H]([C@@H](O)C([O-])=O)C(O)=O', # Tartrate monoanion
        'O[C@@H]([C@H](O)C(O)=O)C(O)=O',  # D-tartaric acid 
        'O[C@@H]([C@H](O)C([O-])=O)C([O-])=O', # D-tartrate dianion
        'O[C@@H]([C@H](O)C([O-])=O)C(O)=O', # D-tartrate monoanion
        'OC(C(O)C([O-])=O)C([O-])=O', # Non-stereochemistry specified
        'OC(C(O)C(O)=O)C(O)=O'
    ]

    # Convert patterns to mol objects
    pattern_mols = [Chem.MolFromSmiles(p) for p in tartrate_patterns]

    tartrate_found = False
    cation_fragments = []
    
    for frag in fragments:
        frag_is_tartrate = False
        for pattern in pattern_mols:
            if frag.HasSubstructMatch(pattern):
                tartrate_found = True
                frag_is_tartrate = True
                break
                
        if not frag_is_tartrate:
            # Check for charged species
            for atom in frag.GetAtoms():
                if atom.GetFormalCharge() > 0 or atom.GetSymbol() == 'N':
                    cation_fragments.append(Chem.MolToSmiles(frag))
                    break
            
    if not tartrate_found:
        return False, "No tartrate anion found"
        
    if not cation_fragments:
        # Look for basic nitrogen atoms that could be protonated
        for frag in fragments:
            if frag_is_tartrate:
                continue
            for atom in frag.GetAtoms():
                if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() > 0:
                    cation_fragments.append(Chem.MolToSmiles(frag))
                    break
                    
    if not cation_fragments:
        return False, "No cation found"
        
    return True, f"Tartrate salt found with cation fragments: {', '.join(cation_fragments)}"
# Pr=1.0
# Recall=1.0