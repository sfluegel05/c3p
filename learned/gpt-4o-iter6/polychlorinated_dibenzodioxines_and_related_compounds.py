"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Classifies a molecule as a polychlorinated dibenzodioxin or structurally related compound based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sufficient number of chlorine or bromine atoms
    halogen_count = sum(atom.GetSymbol() in ['Cl', 'Br'] for atom in mol.GetAtoms())
    if halogen_count < 4:  # Typically related compounds are highly chlorinated
        return False, "Insufficient halogen atoms for a likely related compound."

    # SMARTS pattern for polychlorinated biphenyls (PCBs)
    biphenyl_pattern = Chem.MolFromSmarts('c1cc(Cl)c(Cl)cc1-c2c(Cl)cc(Cl)cc2')
    
    # SMARTS pattern for polychlorinated dibenzodioxins (PCDDs)
    dioxin_patterns = [
        Chem.MolFromSmarts('c1cc2Oc3cc(Cl)cc(Cl)c3Oc2cc1'),
        Chem.MolFromSmarts('c1cc2Oc3ccc(cc3Oc2cc1)Cl')
    ]
    
    # SMARTS pattern for polychlorinated dibenzofurans (PCDFs)
    furan_patterns = [
        Chem.MolFromSmarts('c1cc2oc3cc(Cl)c(Cl)c(Cl)c3c2cc1'),
        Chem.MolFromSmarts('c1cc2oc3ccc(cc3)c2cc1')
    ]
    
    # Check for presence of these patterns in the molecule
    if (any(mol.HasSubstructMatch(pattern) for pattern in dioxin_patterns + furan_patterns + [biphenyl_pattern])):
        return True, "Matches polychlorinated dibenzodioxin or related compound pattern."

    return False, "Does not match polychlorinated dibenzodioxin or related compound pattern."