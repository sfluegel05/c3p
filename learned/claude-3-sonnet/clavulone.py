"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone
Definition: A class of esterified prostanoids obtained from marine corals.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core patterns characteristic of clavulones
    # More specific cyclopentenone core patterns with proper substitution
    core_pattern1 = Chem.MolFromSmarts("[Cl,Br,I,O][C]1[C](=O)[C]=[C][C]1([O,C])[C]") # Halogenated/oxygenated core
    core_pattern2 = Chem.MolFromSmarts("[C]1=C[C]([O,C])[C](=O)[C]1=C") # Conjugated core
    
    # Characteristic side chain patterns
    chain_pattern1 = Chem.MolFromSmarts("CC=CCC") # Alkenyl chain
    chain_pattern2 = Chem.MolFromSmarts("C=CC=C[CH]CC(=O)O[CH3]") # Conjugated ester chain
    
    # Acetoxy and ester groups
    acetoxy = Chem.MolFromSmarts("OC(=O)C")
    ester = Chem.MolFromSmarts("C(=O)OC")
    
    # Check for characteristic core structure
    if not (mol.HasSubstructMatch(core_pattern1) or mol.HasSubstructMatch(core_pattern2)):
        return False, "Missing characteristic clavulone core structure"
    
    # Must have at least one characteristic side chain
    if not (mol.HasSubstructMatch(chain_pattern1) or mol.HasSubstructMatch(chain_pattern2)):
        return False, "Missing characteristic side chain"
    
    # Must have oxygen-containing groups
    acetoxy_count = len(mol.GetSubstructMatches(acetoxy))
    ester_count = len(mol.GetSubstructMatches(ester))
    if acetoxy_count + ester_count < 1:
        return False, "Missing required oxygen-containing groups"
    
    # Check for presence of halogens or oxygen substituents
    hal_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9,17,35,53])
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if hal_count == 0 and o_count < 3:
        return False, "Insufficient characteristic substituents"
    
    # Count carbons to ensure molecule is in right size range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (15 <= c_count <= 30):
        return False, "Carbon count outside typical range for clavulones"
    
    # Additional check for ring count to avoid complex polycyclic compounds
    ring_count = len(Chem.GetSymmSSSR(mol))
    if ring_count > 2:
        return False, "Too many rings for clavulone structure"

    return True, "Contains characteristic clavulone structural features"