from rdkit import Chem
from rdkit.Chem import AllChem
import re

def is_maleate_salt(smiles: str):
    """
    Determines if a molecule is a maleate salt.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a maleate salt, False otherwise
        str: Reason for classification
    """
    try:
        # Check for maleate pattern
        maleate_pattern = r'OC\(=O\)\\C=C/C\(O\)=O|OC\(=O\)/C=C\\C\(O\)=O'
        
        # Check for ionic maleate pattern
        ionic_pattern = r'\[O-\]C\(=O\)\\C=C/C\(\[O-\]\)=O|\[O-\]C\(=O\)/C=C\\C\(\[O-\]\)=O'
        
        if not (re.search(maleate_pattern, smiles) or re.search(ionic_pattern, smiles)):
            return False, "No maleate ion found"

        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Split into fragments
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        
        if len(fragments) < 2:
            return False, "No salt present - single molecule"

        # Look for maleate fragment
        maleate_count = 0
        for frag in fragments:
            frag_smiles = Chem.MolToSmiles(frag)
            if (re.search(maleate_pattern, frag_smiles) or 
                re.search(ionic_pattern, frag_smiles)):
                maleate_count += 1
                
        if maleate_count == 0:
            return False, "No maleate fragment found"
        elif maleate_count == 1:
            return True, "Contains maleate salt"
        else:
            return True, f"Contains {maleate_count} maleate ions (dimaleate/polymaleate salt)"

    except Exception as e:
        return None, f"Error analyzing molecule: {str(e)}"
# Pr=None
# Recall=None