"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for phosphocholine group (-OP(=O)([O-])OCC[N+](C)(C)C)
    phosphocholine = Chem.MolFromSmarts('[OP](=O)([O-])OCC[N+](C)(C)C')
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "Missing phosphocholine group"

    # Check for amide linkage (-NC(=O)-)
    amide = Chem.MolFromSmarts('[NX3][CX3](=[OX1])')
    if not mol.HasSubstructMatch(amide):
        return False, "Missing amide linkage"
    
    # Check for hydroxyl group (characteristic of sphingoid base)
    hydroxyl = Chem.MolFromSmarts('[OH1]')
    if not mol.HasSubstructMatch(hydroxyl):
        return False, "Missing hydroxyl group"

    # Count carbons - sphingomyelins typically have long carbon chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for sphingomyelin"

    # Check molecular weight - sphingomyelins typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for sphingomyelin"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 15:
        return False, "Too few rotatable bonds for sphingomyelin"

    # Check for presence of required atoms
    required_atoms = {'C': 20, 'N': 2, 'O': 4, 'P': 1}
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    
    for element, min_count in required_atoms.items():
        if atom_counts.get(element, 0) < min_count:
            return False, f"Insufficient {element} atoms"

    # Look for sphingoid base pattern (long carbon chain with OH and NH)
    sphingoid_base = Chem.MolFromSmarts('[CH2,CH3]-[CH2]~[CH2]~[CH2]~[CH,CH2]-[CH](-[OH])-[CH](-[NH]-[C])-[CH2]')
    if not mol.HasSubstructMatch(sphingoid_base):
        return False, "Missing sphingoid base pattern"

    return True, "Contains sphingoid base with amide-linked fatty acid and phosphocholine group"