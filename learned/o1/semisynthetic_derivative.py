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

    # Check if molecule is organic (contains carbon atoms)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule is not an organic compound"

    # Calculate molecular weight (no lower limit)
    mol_wt = Descriptors.ExactMolWt(mol)

    # Calculate the number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()

    # Calculate fraction of sp3 carbons
    fraction_csp3 = CalcFractionCSP3(mol)

    # Check for natural product-like features
    # Heuristic: Molecules with at least one ring and moderate-to-high fraction of sp3 carbons
    has_natural_product_features = False
    if num_rings >= 1 and fraction_csp3 >= 0.25:
        has_natural_product_features = True

    # Look for specific natural product-like substructures (e.g., lactone, steroid nucleus)
    # For simplicity, we check for lactone and steroid scaffolds
    lactone_pattern = Chem.MolFromSmarts('O=C1OC([#6])C([#6])C1')
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3CCC4CCCC(C4)C3CCC12')
    if mol.HasSubstructMatch(lactone_pattern) or mol.HasSubstructMatch(steroid_pattern):
        has_natural_product_features = True

    if not has_natural_product_features:
        return False, "Molecule lacks natural product-like features"

    # Check for synthetic modifications
    synthetic_modifications = False

    # Halogenation: Presence of halogen atoms
    if any(atom.GetAtomicNum() in [9, 17, 35, 53, 85] for atom in mol.GetAtoms()):
        synthetic_modifications = True  # Halogenation

    # Acylation: Presence of acyl groups
    acyl_pattern = Chem.MolFromSmarts('C(=O)[#6]')
    if mol.HasSubstructMatch(acyl_pattern):
        synthetic_modifications = True

    # Etherification and Esterification: Presence of ether or ester bonds
    ether_pattern = Chem.MolFromSmarts('C-O-C')
    ester_pattern = Chem.MolFromSmarts('C(=O)O[#6]')
    if mol.HasSubstructMatch(ether_pattern) or mol.HasSubstructMatch(ester_pattern):
        synthetic_modifications = True

    # Methylation and other alkylations: Presence of methyl or alkyl groups attached to heteroatoms
    alkylation = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [7, 8, 16]:  # N, O, S
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and len(neighbor.GetNeighbors()) == 1:
                    alkylation = True  # Alkyl group attached to heteroatom
                    break
    if alkylation:
        synthetic_modifications = True

    # Protecting groups: Common synthetic protecting groups
    boc_pattern = Chem.MolFromSmarts('C(=O)OC(C)(C)C')
    fmoc_pattern = Chem.MolFromSmarts('OCC(=O)O')
    if mol.HasSubstructMatch(boc_pattern) or mol.HasSubstructMatch(fmoc_pattern):
        synthetic_modifications = True

    # Introduction of unnatural functional groups
    azide_pattern = Chem.MolFromSmarts('N=[N+]=[N-]')
    alkyne_pattern = Chem.MolFromSmarts('C#C')
    if mol.HasSubstructMatch(azide_pattern) or mol.HasSubstructMatch(alkyne_pattern):
        synthetic_modifications = True

    if not synthetic_modifications:
        return False, "No synthetic modifications detected"

    return True, "Molecule has natural product-like features with evidence of synthetic modifications"