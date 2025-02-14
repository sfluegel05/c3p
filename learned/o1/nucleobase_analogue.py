"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: CHEBI:25438 nucleobase analogue
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for purine and pyrimidine cores
    purine_core = Chem.MolFromSmarts("c1ncnc2ncnn12")  # Purine core
    pyrimidine_core = Chem.MolFromSmarts("c1cncnc1")    # Pyrimidine core
    
    # Define SMARTS patterns for modified nucleobases
    nucleobase_analogues_patterns = [
        # Modified purines
        Chem.MolFromSmarts("c1ncnc2nc[nH]c12"),       # 2-aminopurine
        Chem.MolFromSmarts("c1ncnc2nc[nH]c(=O)n12"),  # Hypoxanthine analogues
        Chem.MolFromSmarts("c1ncnc2nc[nH]nc12"),      # Azaadenine analogues
        Chem.MolFromSmarts("c1ncnc2nc(=O)[nH]c1n2"),  # 8-oxoadenine analogues
        # Modified pyrimidines
        Chem.MolFromSmarts("O=C1C=CN=C(N1)[O,N,S]"),  # Uracil analogues with substitutions
        Chem.MolFromSmarts("O=C1NC(=O)C=C[N,O,S]1"),  # Thymine analogues with substitutions
        Chem.MolFromSmarts("O=C1NC(=O)C=CN1[O,N,S]"), # Modified uracil/thymine
        Chem.MolFromSmarts("c1nc(N)[nH]cc1[O,N,S]"),  # Cytosine analogues with substitutions
        # Heterocyclic analogues
        Chem.MolFromSmarts("c1n[c,n]c[nH]c1=O"),      # Azauracil analogues
        Chem.MolFromSmarts("c1nc2[nH]nnc2c(=O)[nH]1"),# 8-azaguanine analogues
        # Other analogues
        Chem.MolFromSmarts("c1ncnc2ncnc(=O)c12"),     # Xanthine analogues
        Chem.MolFromSmarts("c1ccn[c,n]c1=O"),         # 5-substituted uracils
        Chem.MolFromSmarts("c1ccn[c,n]c1=O"),         # 5-fluorouracil and analogues
    ]
    
    # Check for purine or pyrimidine core
    purine = mol.HasSubstructMatch(purine_core)
    pyrimidine = mol.HasSubstructMatch(pyrimidine_core)
    
    # Check for matches with nucleobase analogue patterns
    analogue_match = False
    for pattern in nucleobase_analogues_patterns:
        if mol.HasSubstructMatch(pattern):
            analogue_match = True
            break
    
    if purine or pyrimidine or analogue_match:
        return True, "Matches nucleobase analogue pattern"
    
    return False, "Does not match nucleobase analogue patterns"