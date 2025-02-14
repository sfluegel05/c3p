"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: CHEBI:51841 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid is a sesquiterpenoid containing a dihydroagarofuran skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for dihydroagarofuran core
    dihydroagarofuran_core = Chem.MolFromSmarts("[C@H]1[C@]2([C@H](O)[C@@H]([C@H](O)[C@@H]1O)[C@@H]2O)O")
    if not mol.HasSubstructMatch(dihydroagarofuran_core):
        return False, "No dihydroagarofuran core found"

    # Check for ring system and double bond equivalents (DBE) typical of sesquiterpenoids
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 2 or num_rings > 4:
        return False, "Number of rings not typical for sesquiterpenoids"

    dbe = rdMolDescriptors.CalcNumAromaticRings(mol) + rdMolDescriptors.CalcNumRotatableBonds(mol)
    if dbe < 3 or dbe > 7:
        return False, "Double bond equivalents not typical for sesquiterpenoids"

    # Check for common modifications/substituents (e.g., acetylation, benzoylation)
    has_modifications = any(atom.GetSmarts() in ["[#6]-[#6]-[#8]-[#6]", "[#6]-[#6]-[#8]-[#6]-[#6]"] for atom in mol.GetAtoms())

    # Sesquiterpenoid with dihydroagarofuran core, typical ring system and DBE, allowing for common modifications
    return True, "Contains the dihydroagarofuran core and other structural features typical of dihydroagarofuran sesquiterpenoids"