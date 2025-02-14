"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid has a specific bicyclic core with an ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
                         Reason for the classification
    """

    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Core Skeleton Matching using relaxed SMARTS
    # This SMARTS matches the dihydroagarofuran core focusing on the ring connectivity and ether linkage,
    # allowing for substitutions on carbon positions.
    core_pattern = Chem.MolFromSmarts('[C]12[C]([C])[C]([C])[C]([C])([O]1)[C][C]2')
    match = mol.GetSubstructMatch(core_pattern)
    if not match:
        return False, "No dihydroagarofuran core skeleton found."
    
    
    # 3. At least 1 ester group
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
         return False, f"Less than 1 ester groups. Found {len(ester_matches)}"
    
    # 4. Additional checks: At least 4 oxygens, at least 14 carbons
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if oxygen_count < 4:
        return False, f"Too few oxygens. Found {oxygen_count}"
    if carbon_count < 14 :
        return False, f"Too few carbons. Found {carbon_count}"

    # 5. If all checks pass, it's classified as a dihydroagarofuran sesquiterpenoid.
    return True, "Contains a dihydroagarofuran core, at least one ester group, four oxygens and 14 carbons"