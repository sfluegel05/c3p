"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem

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

    # Patterns representing known nucleobase-like cores
    nucleobase_patterns = [
        Chem.MolFromSmarts("c1ncnc2[nH]ncc12"),  # Purine core
        Chem.MolFromSmarts("c1[nH]cc(=O)nc1"),  # Uracil-like core
        Chem.MolFromSmarts("c1c(=O)[nH]cnc1"),  # Cytosine-like core
    ]
    
    # Check for nucleobase-like structures
    for nb_pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(nb_pattern):
            # Further check for modifications (e.g., methylation, halogenation, etc.)
            methyl_pattern = Chem.MolFromSmarts("[CH3]")
            halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
            hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
            keto_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")

            if any(mol.HasSubstructMatch(p) for p in [methyl_pattern, halogen_pattern, hydroxyl_pattern, keto_pattern]):
                return True, "Contains nucleobase-like core with typical modifications"

    # Expanding further to consider known modifications may capture additional nucleobase analogues
    pyrimidine_like_pattern = Chem.MolFromSmarts("c1cnc[nH]c1")
    if mol.HasSubstructMatch(pyrimidine_like_pattern):
        if any(mol.HasSubstructMatch(p) for p in [methyl_pattern, keto_pattern]):
            return True, "Contains pyrimidine-like structure with characteristic modifications"
        
    return False, "Does not match features of known nucleobase analogues"