"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid ring system (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~[#6]4~[#6]~[#6]3~[#6]~[#6]2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid ring system found"

    # Check for hydroxyl group
    hydroxyl = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl):
        return False, "No hydroxyl group found"
    
    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Check carbon count (sterols typically have 27-30 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for a sterol"
    if c_count > 35:
        return False, "Too many carbons for a sterol"

    # Check for characteristic methyl groups at C-10 and C-13
    angular_methyls = Chem.MolFromSmarts("[CH3]~[C]([#6])([#6])[#6]")
    if len(mol.GetSubstructMatches(angular_methyls)) < 1:
        return False, "Missing characteristic methyl groups"

    # Calculate molecular weight (typical sterol MW range)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 600:
        return False, f"Molecular weight {mol_wt:.1f} outside typical sterol range"

    # Check for aliphatic side chain
    side_chain = Chem.MolFromSmarts("[CH2,CH3,CH]~[CH2,CH3,CH]~[CH2,CH3,CH]")
    if not mol.HasSubstructMatch(side_chain):
        return False, "No characteristic side chain found"

    # Count oxygen atoms (sterols typically have 1-3 oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1 or o_count > 6:
        return False, f"Number of oxygen atoms ({o_count}) outside typical range for sterols"

    # Specific check for 3-hydroxy position
    hydroxy_3_pos = Chem.MolFromSmarts("[#6]1~[#6]~[#6]([OH])~[#6]~[#6]~[#6]~1")
    if not mol.HasSubstructMatch(hydroxy_3_pos):
        return False, "No hydroxyl group in characteristic 3-position"

    return True, "Contains steroid ring system with 3-hydroxy group and characteristic sterol features"