from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import rdMolDescriptors

def is_precorrin(smiles: str):
    """
    Determines if a molecule is a precorrin based on structural characteristics.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a precorrin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Key characteristics of precorrins:
    # 1. Contains tetrapyrrole framework
    # 2. Has carboxylic acid groups
    # 3. Contains methyl substituents
    # 4. Has specific molecular weight range
    # 5. Contains nitrogen atoms in ring system
    
    # Check for presence of carboxylic acid groups
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid groups found"
        
    # Count number of carboxylic acid groups (should have multiple)
    num_carboxylic = len(mol.GetSubstructMatches(carboxylic_pattern))
    if num_carboxylic < 4:
        return False, "Insufficient number of carboxylic acid groups"
        
    # Check for pyrrole rings
    pyrrole_pattern = Chem.MolFromSmarts('[nH]1cccc1')
    if not mol.HasSubstructMatch(pyrrole_pattern):
        pyrrole_pattern2 = Chem.MolFromSmarts('n1cccc1') # Check for deprotonated pyrroles too
        if not mol.HasSubstructMatch(pyrrole_pattern2):
            return False, "No pyrrole rings found"
        
    # Count nitrogen atoms (should have 4 or more for tetrapyrrole)
    num_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if num_nitrogens < 4:
        return False, "Insufficient number of nitrogen atoms for tetrapyrrole structure"
        
    # Check molecular weight (typical range for precorrins)
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if not (800 < mol_weight < 1000):
        return False, f"Molecular weight {mol_weight:.1f} outside typical precorrin range"
        
    # Check for methyl groups
    methyl_pattern = Chem.MolFromSmarts('C[CH3]')
    num_methyls = len(mol.GetSubstructMatches(methyl_pattern))
    if num_methyls < 2:
        return False, "Insufficient number of methyl groups"
        
    # Check for ring system size
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 4:
        return False, "Insufficient number of rings for tetrapyrrole structure"
        
    return True, "Molecule contains characteristic precorrin features: tetrapyrrole framework, multiple carboxylic acids, methyl groups"
# Pr=None
# Recall=None