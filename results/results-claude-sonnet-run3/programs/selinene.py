from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_selinene(smiles: str):
    """
    Determines if a molecule is a selinene (sesquiterpene with carbobicyclic skeleton).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a selinene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check molecular formula - should be C15H24 or C15H26O
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if formula not in ["C15H24", "C15H26O"]:
        return False, f"Incorrect molecular formula {formula} - selinenes are C15H24 or C15H26O"
        
    # Check for bicyclic structure
    rings = mol.GetRingInfo()
    if rings.NumRings() != 2:
        return False, f"Not bicyclic - has {rings.NumRings()} rings"
        
    # Check ring sizes - should be 6,6 or 6,7 membered rings
    ring_sizes = sorted([len(ring) for ring in rings.AtomRings()])
    if ring_sizes != [6,6]:
        return False, f"Incorrect ring sizes {ring_sizes} - should be 6,6"
        
    # Check for presence of at least one double bond or hydroxyl group
    double_bond = Chem.MolFromSmarts("C=C")
    hydroxyl = Chem.MolFromSmarts("OH")
    if not (mol.HasSubstructMatch(double_bond) or mol.HasSubstructMatch(hydroxyl)):
        return False, "No double bonds or hydroxyl groups present"
        
    # Check for presence of methyl groups
    pattern = Chem.MolFromSmarts("[CH3]")
    if len(mol.GetSubstructMatches(pattern)) < 2:
        return False, "Less than 2 methyl groups present"
        
    # Check for gem-dimethyl group
    gem_dimethyl = Chem.MolFromSmarts("C([CH3])([CH3])")
    if not mol.HasSubstructMatch(gem_dimethyl):
        return False, "Missing characteristic gem-dimethyl group"

    # Check that rings are fused (share exactly 2 atoms)
    ring_atoms = list(rings.AtomRings())
    shared_atoms = set(ring_atoms[0]).intersection(set(ring_atoms[1]))
    if len(shared_atoms) != 2:
        return False, "Rings are not properly fused"

    # Check for characteristic selinene core structure - decalin with at least one methyl group
    selinene_core = Chem.MolFromSmarts("C1CCC2CCCC2C1")
    if not mol.HasSubstructMatch(selinene_core):
        return False, "Missing characteristic decalin core structure"

    return True, "Matches selinene structure: bicyclic sesquiterpene with correct ring sizes and substitution pattern"
# Pr=None
# Recall=0.0