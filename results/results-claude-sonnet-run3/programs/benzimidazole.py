from rdkit import Chem
from rdkit.Chem import AllChem

def is_benzimidazole(smiles: str):
    """
    Determines if a molecule is a benzimidazole or contains a benzimidazole core.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is/contains benzimidazole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # SMARTS patterns for benzimidazole core variants
    patterns = [
        # 1H-benzimidazole
        "c1ccc2[nH]cnc2c1",
        # 2H-benzimidazole 
        "C1N=c2ccccc2=N1",
        # 3aH-benzimidazole
        "C1=CC2N=CN=C2C=C1",
        # 4H-benzimidazole
        "C1C=CC=C2N=CN=C12"
    ]

    for pattern in patterns:
        patt = Chem.MolFromSmiles(pattern)
        if mol.HasSubstructMatch(patt):
            return True, "Contains benzimidazole core structure"

    # Check for basic ring structure
    rings = mol.GetRingInfo()
    
    # Need at least 2 rings
    if rings.NumRings() < 2:
        return False, "Does not contain required ring system"

    # Look for fused 6 and 5 membered rings
    has_6ring = False
    has_5ring = False
    
    for ring in rings.AtomRings():
        if len(ring) == 6:
            # Check if 6-membered ring atoms are all carbon
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetSymbol() == 'C' for atom in atoms):
                has_6ring = True
        if len(ring) == 5:
            # Check if 5-membered ring has 2 nitrogens
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            n_count = sum(1 for atom in atoms if atom.GetSymbol() == 'N')
            if n_count == 2:
                has_5ring = True
            
    if not (has_5ring and has_6ring):
        return False, "Does not contain required benzene and imidazole ring system"

    return False, "Does not match benzimidazole structure pattern"
# Pr=0.2727272727272727
# Recall=1.0