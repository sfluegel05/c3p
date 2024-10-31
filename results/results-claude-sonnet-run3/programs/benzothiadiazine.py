from rdkit import Chem
from rdkit.Chem import AllChem

def is_benzothiadiazine(smiles: str):
    """
    Determines if a molecule is a benzothiadiazine.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a benzothiadiazine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for presence of required atoms
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if not ('S' in atoms and atoms.count('N') >= 2):
        return False, "Missing required S and/or N atoms"

    # Find benzene ring
    benzene = mol.GetSubstructMatches(Chem.MolFromSmarts('c1ccccc1'))
    if not benzene:
        return False, "No benzene ring found"

    # Check for SO2 group
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('S(=O)(=O)')):
        return False, "Missing SO2 group"

    # SMARTS patterns for different benzothiadiazine core variations
    patterns = [
        # Basic benzothiadiazine core with SO2 group
        '[#6]1:[#6]:[#6]:[#6]:[#6]2:[#6]:1[#16](=[#8])(=[#8])[#7]~[#6]~[#7]~2',
        # Variation with single bonds in N-containing ring
        '[#6]1:[#6]:[#6]:[#6]:[#6]2:[#6]:1[#16](=[#8])(=[#8])[#7]-[#6]-[#7]-2',
        # Variation with double bond to one N
        '[#6]1:[#6]:[#6]:[#6]:[#6]2:[#6]:1[#16](=[#8])(=[#8])[#7]=[#6]-[#7]-2',
        # Variation with double bond between Ns
        '[#6]1:[#6]:[#6]:[#6]:[#6]2:[#6]:1[#16](=[#8])(=[#8])[#7]-[#6]=[#7]-2'
    ]

    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            return True, "Contains benzothiadiazine core structure"

    return False, "Missing benzothiadiazine core structure"
# Pr=1.0
# Recall=0.9565217391304348