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
    
    # Define valid core structure SMARTS patterns
    # Dibenzodioxin: two benzene rings connected by two oxygen atoms
    dioxin_core = Chem.MolFromSmarts("O1c2ccccc2Oc3ccccc13")
    # Dibenzofuran: two benzene rings connected by one oxygen atom
    furan_core = Chem.MolFromSmarts("O1c2ccccc2c3ccccc13")
    # Biphenyl core (for PCBs/PBBs)
    biphenyl_core = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    
    # Verify SMARTS compilation
    core_structures = [dioxin_core, furan_core, biphenyl_core]
    if any(core is None for core in core_structures):
        return None, None  # Handle invalid SMARTS patterns
    
    # Check for any core structure match
    core_found = any(mol.HasSubstructMatch(core) for core in core_structures)
    if not core_found:
        return False, "No dioxin, furan, or biphenyl core found"
    
    # Count total chlorine and bromine atoms
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35])
    
    if halogen_count < 2:
        return False, f"Only {halogen_count} halogen atoms (needs ≥2)"
    
    return True, f"Contains core structure with {halogen_count} halogen atoms (Cl/Br)"