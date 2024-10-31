from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sapogenin(smiles: str):
    """
    Determines if a molecule is a sapogenin based on structural characteristics.
    Sapogenins are polycyclic compounds that are the aglycon part of saponins,
    typically with steroid or triterpenoid skeletons.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for basic steroid/triterpenoid core patterns
    steroid_pattern = Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1')
    triterpenoid_pattern = Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1')
    
    if not (mol.HasSubstructMatch(steroid_pattern) or mol.HasSubstructMatch(triterpenoid_pattern)):
        return False, "Missing steroid/triterpenoid core structure"

    # Count rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 4:
        return False, "Too few rings for a sapogenin"

    # Check for oxygen-containing groups
    hydroxyl_pattern = Chem.MolFromSmarts('[OH]')
    ether_pattern = Chem.MolFromSmarts('[OR]')
    spiro_oxygen_pattern = Chem.MolFromSmarts('O1CC[CH2,CH]C1')
    
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)
    has_spiro_oxygen = mol.HasSubstructMatch(spiro_oxygen_pattern)

    if not (has_hydroxyl or has_ether or has_spiro_oxygen):
        return False, "Missing required oxygen-containing groups"

    # Count carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons < 20:
        return False, "Too few carbons for a sapogenin"

    features = []
    if has_hydroxyl:
        features.append("hydroxyl group")
    if has_spiro_oxygen:
        features.append("spiro oxygen ring")
    if has_ether:
        features.append("ether linkage")
    features.append(f"{num_rings} rings")
    
    return True, f"Sapogenin structure confirmed: {', '.join(features)}"
# Pr=None
# Recall=0.0