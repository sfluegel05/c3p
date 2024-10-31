from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy steroid.
    A 3-hydroxy steroid is defined as any steroid carrying a hydroxy group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # First check for steroid core
        steroid_core = Chem.MolFromSmarts('[C]1[C][C]2[C]3[C][C][C]4[C@@H]([C@@H]3[C@H]([C@H]([C]2)1)[C,c])[C,c][C,c][C,c]4')
        if not mol.HasSubstructMatch(steroid_core):
            return False, "No steroid core found"

        # SMARTS patterns for 3-hydroxy group in steroids
        patterns = [
            # 3-hydroxy with tetrahedral carbon
            '[OH1][C@H,@@H]1[CH2][CH2][C,c]2[C,c][C,c][C,c]',
            # 3-hydroxy with aromatic A-ring
            '[OH1]c1cc[c,C]2[C,c][C,c][C,c]',
            # More general 3-hydroxy pattern
            '[OH1][C,c]1[C,c][C,c][C,c]2[C,c][C,c][C,c]'
        ]

        for pattern in patterns:
            patt = Chem.MolFromSmarts(pattern)
            if patt and mol.HasSubstructMatch(patt):
                # Additional check to confirm position 3
                matches = mol.GetSubstructMatches(patt)
                for match in matches:
                    # The OH group should be attached to the third carbon of the first ring
                    if any(idx in match[:2] for idx in mol.GetSubstructMatch(steroid_core)[:3]):
                        return True, "Found hydroxyl group at position 3 of steroid core"

        return False, "No hydroxyl group at position 3"

    except Exception as e:
        return None, f"Error analyzing molecule: {str(e)}"
# Pr=None
# Recall=0.0