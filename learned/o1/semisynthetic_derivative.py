"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticRings, CalcFractionCSP3

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    A semisynthetic derivative is an organic molecule derived from a natural product by partial chemical synthesis.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is organic
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule is not an organic compound"

    # Calculate molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low to be a typical semisynthetic derivative"

    # Calculate the number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return False, "No ring structures found, unlikely to be a natural product derivative"

    # Calculate fraction of sp3 carbons
    fraction_csp3 = CalcFractionCSP3(mol)
    if fraction_csp3 < 0.25:
        return False, "Low fraction of sp3 carbons, less likely to be derived from natural products"

    # Check for presence of natural product-like scaffold
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smiles = Chem.MolToSmiles(scaffold)
    # Heuristic: If scaffold is complex (contains multiple rings)
    scaffold_ring_info = scaffold.GetRingInfo()
    scaffold_num_rings = scaffold_ring_info.NumRings()
    if scaffold_num_rings < 2:
        return False, "Scaffold is too simple, unlikely to be from a natural product"

    # Look for synthetic modifications
    # Common modifications: acetylation, methylation, halogenation, esterification
    synthetic_modifications = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [9, 17, 35, 53]:
            synthetic_modifications = True  # Halogenation
            break
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    if mol.HasSubstructMatch(ester_pattern):
        synthetic_modifications = True  # Esterification
    acetyl_pattern = Chem.MolFromSmarts('C(=O)C')
    if mol.HasSubstructMatch(acetyl_pattern):
        synthetic_modifications = True  # Acetylation
    methylation = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and len(atom.GetNeighbors()) == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetAtomicNum() == 6:
                methylation = True  # Methyl group attached to carbon
                break
    if methylation:
        synthetic_modifications = True

    if not synthetic_modifications:
        return False, "No synthetic modifications detected"

    return True, "Molecule has a natural product-like scaffold with synthetic modifications"