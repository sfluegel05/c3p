"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids with specific structural characteristics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count Rings (must be 4 or more)
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings < 4:
        return False, "Less than 4 rings"
    
    # Molecular Weight check - cucurbitacins are triterpenoids, usually > 400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 700:
      return False, f"Molecular weight {mol_wt} is out of the typical range (400-700 Da)"

    # Count Carbons (must be about 30, but we will check the mol_wt instead )
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27:
        return False, f"Too few carbons {c_count}, expected at least 27"

    # Count Oxygens (must be at least 4)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
         return False, f"Oxygen count {o_count} is too low (must be at least 3)"
    
    # Check for tetracyclic core using a more specific SMARTS pattern
    # This pattern looks for a fused 4-ring system
    core_pattern = Chem.MolFromSmarts("C1[C]2[C]3[C]4[C]1[C]5[C]2[C]3[C]45") 
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Tetracyclic core not found"

    # Check for at least 3 oxygens or carbonyls on the ring system (more than just alcohol groups)
    # Check for common side chain patterns
    # This includes alcohols and carbonyls attached to the core rings
    oxygen_pattern = Chem.MolFromSmarts("[CX4]([OX2])[CX4]") # hydroxyls on carbons attached to rings
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]") # carbonyl groups
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)

    ring_oxygen_pattern = Chem.MolFromSmarts("[O;R]") #oxygens within the rings
    ring_oxygen_matches = mol.GetSubstructMatches(ring_oxygen_pattern)

    if (len(oxygen_matches) + len(carbonyl_matches) + len(ring_oxygen_matches)) < 3:
       return False, "Not enough oxygens on the core rings"

     # Check for possible isopropyl group attached to a ring carbon
    isopropyl_pattern = Chem.MolFromSmarts("[C;R][CX4](C)(C)")
    isopropyl_matches = mol.GetSubstructMatches(isopropyl_pattern)
    if len(isopropyl_matches) < 1 : # many cucurbitacins have the isopropyl, some don't.
          pass # don't fail, as isopropyl isn't mandatory, and may be substituted

    # Check for a double bond outside the ring structure
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 1:
      pass # double bonds are common, but not necessary.
          
    # If all checks passed, we consider it a cucurbitacin
    return True, "Likely a cucurbitacin based on structure features"