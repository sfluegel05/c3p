"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
Definition: A flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic molecular properties
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 30:  # Minimum size for dimeric structure
        return False, "Too small to be a proanthocyanidin"

    # Count key atoms
    num_o = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_o < 8:  # Minimum oxygens for dimeric structure
        return False, "Insufficient oxygen atoms"

    # Flavan-3-ol unit patterns
    flavan_patterns = [
        # Basic flavan-3-ol core (both 2,3-stereochemistries)
        "[#6]1-[#6@H]-[#6@H]-c2c(O)cc(O)cc2O[#6H]-1",
        "[#6]1-[#6@@H]-[#6@H]-c2c(O)cc(O)cc2O[#6H]-1",
        # Gallocatechin variants
        "[#6]1-[#6@H]-[#6@H]-c2c(O)cc(O)cc2O[#6H]-1-c1cc(O)c(O)c(O)c1",
        "[#6]1-[#6@@H]-[#6@H]-c2c(O)cc(O)cc2O[#6H]-1-c1cc(O)c(O)c(O)c1"
    ]
    
    total_flavan_matches = 0
    for pattern in flavan_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt:
            matches = len(mol.GetSubstructMatches(patt))
            total_flavan_matches += matches

    if total_flavan_matches < 2:
        return False, "Insufficient flavan-3-ol units"

    # Proanthocyanidin linkage patterns
    linkage_patterns = [
        # 4→8 linkage (both stereochemistries)
        "[#6]1-[#6@H]-[#6@H]-c2c(-[#6@H]1)cc(O)c1c2O[#6]-[#6]-[#6]c2cc(O)cc(O)c21",
        "[#6]1-[#6@@H]-[#6@H]-c2c(-[#6@H]1)cc(O)c1c2O[#6]-[#6]-[#6]c2cc(O)cc(O)c21",
        # 4→6 linkage
        "[#6]1-[#6@H]-[#6@H]-c2c(-[#6@H]1)cc(O)c(-[#6]1-[#6]-[#6]-c3c(-[#6]1)cc(O)cc3O)c2O",
        # A-type (2→O→7, 4→8)
        "[#6]1O[#6]2Oc3cc(O)cc(O)c3[#6@H][#6@H]2[#6][#6]1"
    ]
    
    has_linkage = False
    for pattern in linkage_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            has_linkage = True
            break

    if not has_linkage:
        return False, "Missing characteristic proanthocyanidin linkages"

    # Required substructure patterns
    required_patterns = [
        # A-ring with meta-hydroxylation
        "Oc1cc(O)cc(C2)c1O[C@@H]([C@H]2O)",
        # B-ring patterns
        "c1cc(O)c(O)cc1",  # catechin type
        "c1cc(O)c(O)c(O)c1",  # gallocatechin type
    ]
    
    missing_patterns = 0
    for pattern in required_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and not mol.HasSubstructMatch(patt):
            missing_patterns += 1

    if missing_patterns > 1:  # Allow some flexibility
        return False, "Missing required structural features"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 550:  # Minimum weight for dimeric structure
        return False, "Molecular weight too low for proanthocyanidin"

    # Count aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    if aromatic_rings < 4:  # Minimum for dimeric structure
        return False, "Insufficient aromatic rings"

    # Final check for oligomeric nature
    if total_flavan_matches >= 2 and has_linkage:
        return True, "Contains multiple flavan-3-ol units with characteristic proanthocyanidin linkages"
    
    return False, "Does not meet structural requirements for proanthocyanidin"