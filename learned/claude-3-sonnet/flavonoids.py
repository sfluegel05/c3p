"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: CHEBI:44102 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a phenyl-substituted 1-phenylpropane skeleton
    with a C15 or C16 skeleton, or a structure condensed with a C6-C3 lignan precursor.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define common flavonoid substructure patterns
    flavonoid_patterns = [
        Chem.MolFromSmarts("[c1c(O)cc(O)cc1]"),  # Catechol ring
        Chem.MolFromSmarts("[c1c(O)ccc(O)c1]"),  # Resorcinol ring
        Chem.MolFromSmarts("[C1=C(O)C(=O)C=C(O)C1]"),  # Chromen-4-one
        Chem.MolFromSmarts("[C1=C(O)C(=O)C2=C(O)C=CC=C2O1]"),  # Flavone
        Chem.MolFromSmarts("[C1=C(O)C(=O)C2=C(O)C=CC(O)=C2O1]"),  # Flavonol
        Chem.MolFromSmarts("[C1=CC(=O)C2=C(O)C=C(O)C=C2O1]"),  # Isoflavone
        Chem.MolFromSmarts("[C1=C(O)C(=O)C2=C(O)C=C(O)C=C2O1]"),  # Anthocyanidin
        Chem.MolFromSmarts("[C1=C(O)C(=O)C2=C(O)C=CC=C2O1]"),  # Flavanone
        Chem.MolFromSmarts("[C1=C(O)C(=O)C2=C(O)C=CC(O)=C2O1]"),  # Dihydroflavonol
        Chem.MolFromSmarts("[C1=C(O)C(=O)C2=C(O)C=C(O)C=C2O1]"),  # Flavonol
    ]
    
    # Check for flavonoid substructure patterns
    has_flavonoid_pattern = any(mol.HasSubstructMatch(pattern) for pattern in flavonoid_patterns)
    if not has_flavonoid_pattern:
        return False, "No common flavonoid substructure found"
    
    # Check for phenyl-substituted 1-phenylpropane skeleton or C6-C3 lignan precursor
    phenylpropane_pattern = Chem.MolFromSmarts("[c1ccccc1]CCC[c2]ccccc2")
    lignan_pattern = Chem.MolFromSmarts("[c1ccccc1]CC[c2]ccccc2")
    has_skeleton = mol.HasSubstructMatch(phenylpropane_pattern) or mol.HasSubstructMatch(lignan_pattern)
    if not has_skeleton:
        return False, "No flavonoid skeleton found"
    
    # Check for C15 or C16 skeleton
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count not in [15, 16]:
        return False, "Carbon skeleton size not C15 or C16"
    
    # Check for oxygen count (typically 3-7 oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3 or o_count > 7:
        return False, "Oxygen count outside typical range for flavonoids"
    
    # Check for aromatic rings (typically 2-4 rings)
    ring_info = mol.GetRingInfo()
    n_aromatic_rings = sum(ring_info.IsAromaticRing(i) for i in range(ring_info.NumRings()))
    if n_aromatic_rings < 2 or n_aromatic_rings > 4:
        return False, "Number of aromatic rings outside typical range for flavonoids"
    
    # Check for molecular weight (typically 200-800 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, "Molecular weight outside typical range for flavonoids"
    
    return True, "Molecule possesses structural features consistent with flavonoids"