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

    # 2. Core Skeleton Matching using SMARTS
    # This SMARTS matches the dihydroagarofuran core, with * as variable carbon positions.
    # Includes the methyl group on the core, and focuses on the ring connectivity. Ignores stereochemistry for now.
    core_pattern = Chem.MolFromSmarts('[C]12[C]([C][O]1)[C][C][C]([C])([C]2[C])')
    match = mol.GetSubstructMatch(core_pattern)
    if not match:
        return False, "No dihydroagarofuran core skeleton found."
    
    # 3. Stereochemistry check for core atoms.
    # Get the indices for the core atoms
    c1, c2, c3, c4, c5, c6, c7, c8  = match
    core_atoms = [c1,c2,c3,c4,c5,c6,c7,c8]

    #Check stereochemistry at C1 and C6 (bridgehead carbons) and C2
    c1_atom = mol.GetAtomWithIdx(c1)
    c6_atom = mol.GetAtomWithIdx(c6)
    c2_atom = mol.GetAtomWithIdx(c2)

    if not (c1_atom.HasProp('_CIPCode') and c6_atom.HasProp('_CIPCode') and c2_atom.HasProp('_CIPCode')):
        AllChem.AssignStereochemistry(mol, cleanIt=True, force=True)

    if c1_atom.HasProp('_CIPCode') and c6_atom.HasProp('_CIPCode') and c2_atom.HasProp('_CIPCode') :
        if (c1_atom.GetProp('_CIPCode') != 'R' ) or (c6_atom.GetProp('_CIPCode') != 'R') or (c2_atom.GetProp('_CIPCode') != 'S'):
            return False, "Incorrect stereochemistry at core atoms."
    else:
        return False, "Could not determine stereochemistry of the core"


    # 4. At least 2 ester groups
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
         return False, f"Less than 2 ester groups. Found {len(ester_matches)}"
    
    # 5. Additional checks: At least 4 oxygens, at least 15 carbons
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if oxygen_count < 4:
        return False, f"Too few oxygens. Found {oxygen_count}"
    if carbon_count < 15 :
        return False, f"Too few carbons. Found {carbon_count}"

    # 6. If all checks pass, it's classified as a dihydroagarofuran sesquiterpenoid.
    return True, "Contains a dihydroagarofuran core with correct stereochemistry, at least two ester groups, four oxygens and 15 carbons"