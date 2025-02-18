"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Polychlorinated dibenzodioxins and related compounds (persistent organic pollutants)
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxins and related compounds.
    These include polychlorinated dibenzodioxins, dibenzofurans, biphenyls (PCBs), and polybrominated biphenyls (PBBs).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule matches the class, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define core structure SMARTS patterns
    dioxin_core = Chem.MolFromSmarts("c12c3ccccc3Oc4ccccc4c1")  # Dibenzodioxin structure
    furan_core = Chem.MolFromSmarts("c12c3ccccc3Oc4cccc1c4")    # Dibenzofuran structure
    biphenyl_core = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")     # Biphenyl structure
    
    # Check for any of the core structures
    core_found = (mol.HasSubstructMatch(dioxin_core) or 
                  mol.HasSubstructMatch(furan_core) or 
                  mol.HasSubstructMatch(biphenyl_core))
    
    if not core_found:
        return False, "No dioxin, furan, or biphenyl core found"
    
    # Count total chlorine and bromine atoms
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35])
    
    if halogen_count < 2:
        return False, f"Only {halogen_count} halogen atoms (needs ≥2)"
    
    return True, "Contains core structure with ≥2 halogen atoms (Cl/Br)"