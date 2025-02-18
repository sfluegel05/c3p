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
    These are organochlorine/bromine compounds with dioxin, dibenzofuran, or biphenyl cores,
    with halogens attached to aromatic carbons.
    
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

    # Define core structure patterns with accurate SMARTS
    dioxin_core = Chem.MolFromSmarts('O1c2ccccc2Oc3ccccc13')  # Dibenzodioxin scaffold
    furan_core = Chem.MolFromSmarts('c1ccc2c(c1)oc3ccccc23')  # Dibenzofuran scaffold
    biphenyl_core = Chem.MolFromSmarts('c1ccccc1-c2ccccc2')   # Biphenyl scaffold

    # Check for core structure presence
    core_found = any(mol.HasSubstructMatch(core) for core in [dioxin_core, furan_core, biphenyl_core])
    if not core_found:
        return False, "No dioxin/dibenzofuran/biphenyl core found"

    # Check halogens are attached to aromatic carbons
    halogen_attached = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [17, 35]:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic() and neighbor.GetAtomicNum() == 6:
                    halogen_attached = True
                    break
            if halogen_attached:
                break
    if not halogen_attached:
        return False, "No halogens attached to aromatic carbons"

    # Verify minimum aromatic rings based on core type
    if rdMolDescriptors.CalcNumAromaticRings(mol) < 2:
        return False, "Insufficient aromatic rings for core structure"

    # Additional check for at least 3 halogen atoms (common in POPs)
    if halogen_count < 3:
        return False, f"Insufficient halogens ({halogen_count} found, minimum 3)"

    return True, "Halogenated dioxin/dibenzofuran/biphenyl core with aromatic-bound halogens"