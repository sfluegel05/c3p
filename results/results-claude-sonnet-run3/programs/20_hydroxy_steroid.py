from rdkit import Chem
from rdkit.Chem import AllChem
import re

def is_20_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 20-hydroxy steroid (steroid with a hydroxy group at position 20).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 20-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Steroid core pattern (four fused rings)
        steroid_core = Chem.MolFromSmarts(
            "[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~2~[#6]~1"
        )
        
        if not mol.HasSubstructMatch(steroid_core):
            return False, "No steroid core found"

        # Pattern for 20-hydroxy group
        # Looking for a carbon with OH group connected to the D ring of steroid
        hydroxy_pattern = Chem.MolFromSmarts(
            "[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~2~[#6]~1-[#6](-[OH])"
        )

        # Alternative patterns for 20-hydroxy group
        alt_patterns = [
            "[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~2~[#6]~1[#6](-[OH])",
            "[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~2~[#6]~1-[#6](-[OH])-[#6]",
            "[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~2~[#6]~1[#6](-[OH])-[#6]"
        ]

        if mol.HasSubstructMatch(hydroxy_pattern):
            return True, "Found steroid core with hydroxy group at position 20"

        for pattern in alt_patterns:
            pattern_mol = Chem.MolFromSmarts(pattern)
            if mol.HasSubstructMatch(pattern_mol):
                return True, "Found steroid core with hydroxy group at position 20"

        return False, "No hydroxy group found at position 20"

    except Exception as e:
        return None, f"Error analyzing molecule: {str(e)}"
# Pr=None
# Recall=0.0