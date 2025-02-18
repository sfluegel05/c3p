"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: polychlorinated dibenzodioxines and related compounds (persistent organic pollutants)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to polychlorinated dibenzodioxines and related compounds.
    These are organochlorine/bromine compounds with dioxin, dibenzofuran, or biphenyl cores.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if matches criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Must contain at least one chlorine or bromine
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35])
    if halogen_count == 0:
        return False, "No chlorine/bromine atoms found"

    # Core structure patterns (allowing substitutions)
    dioxin_core = Chem.MolFromSmarts('[O]1c2cccc([O]c3c1cccc3)c2')  # Dibenzodioxin scaffold
    furan_core = Chem.MolFromSmarts('[O]1c2cccc(oc3c1cccc3)c2')    # Dibenzofuran scaffold
    biphenyl_core = Chem.MolFromSmarts('[c]1[c][c][c][c][c]1-[c]2[c][c][c][c][c]2')  # Biphenyl scaffold

    # Check for core structure matches
    core_found = False
    for core in [dioxin_core, furan_core, biphenyl_core]:
        if mol.HasSubstructMatch(core):
            core_found = True
            break
            
    if not core_found:
        return False, "No dioxin/dibenzofuran/biphenyl core found"

    # Check for halogen attachment to aromatic systems
    aromatic_halogens = sum(1 for atom in mol.GetAtoms() 
                          if atom.GetAtomicNum() in [17,35] 
                          and atom.GetIsAromatic())
    if aromatic_halogens == 0:
        return False, "Halogens not attached to aromatic systems"

    # Additional check for proper substitution pattern
    if rdMolDescriptors.CalcNumAromaticRings(mol) < 2:
        return False, "Insufficient aromatic rings for core structure"

    return True, "Contains halogenated dioxin/dibenzofuran/biphenyl core structure"