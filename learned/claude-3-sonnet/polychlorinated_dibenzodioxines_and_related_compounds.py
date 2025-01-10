"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxine or related compound.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count halogens (Cl, Br)
    num_cl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Cl]')))
    num_br = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[Br]')))
    total_halogens = num_cl + num_br
    
    if total_halogens == 0:
        return False, "No chlorine or bromine atoms found"

    # Check for biphenyl core structure
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c1ccccc1')
    
    # Check for dibenzodioxin core structure (two oxygen bridges)
    dioxin_pattern = Chem.MolFromSmarts('c1ccccc1Oc1ccccc1O')
    
    # Check for dibenzofuran core structure (one oxygen bridge)
    furan_pattern = Chem.MolFromSmarts('c1ccccc1oc1ccccc1')
    
    is_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
    is_dioxin = mol.HasSubstructMatch(dioxin_pattern)
    is_furan = mol.HasSubstructMatch(furan_pattern)
    
    # Count aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    # Classify based on structural features
    if is_dioxin and total_halogens >= 1:
        return True, f"Polychlorinated dibenzodioxin with {total_halogens} halogen(s)"
    
    elif is_furan and total_halogens >= 1:
        return True, f"Polychlorinated dibenzofuran with {total_halogens} halogen(s)"
    
    elif is_biphenyl and total_halogens >= 1:
        halogen_type = "chlorinated" if num_cl > 0 else ""
        halogen_type += " and brominated" if num_br > 0 else ""
        return True, f"Poly{halogen_type} biphenyl with {total_halogens} halogen(s)"
    
    # Special case for complex structures like formicamycins that contain these motifs
    elif aromatic_rings >= 2 and total_halogens >= 1:
        # Check if molecule contains both chlorine and oxygen
        has_oxygen = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
        if has_oxygen:
            return True, "Complex polychlorinated aromatic compound with oxygen-containing linkages"
    
    return False, "Does not match structural patterns for PCDDs, PCDFs, or PCBs/PBBs"