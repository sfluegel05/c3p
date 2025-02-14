"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: CHEBI:36326 lipopolysaccharide

Lipopolysaccharides (LPS) are complex natural compounds consisting of a trisaccharide repeating unit 
(two heptose units and octulosonic acid) with oligosaccharide side chains and 3-hydroxytetradecanoic 
acid units. They are a major constituent of the cell walls of Gram-negative bacteria.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.Fingerprints import FingerprintMols

def is_lipopolysaccharide(smiles):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # LPS typically have high molecular weights (>1000 Da)
    if mol_wt < 1000:
        return False, "Molecular weight too low for lipopolysaccharide"
    
    # Calculate oxygen-to-carbon ratio
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_to_c_ratio = o_count / c_count if c_count > 0 else 0
    
    # LPS typically have high oxygen-to-carbon ratio (>0.5)
    if o_to_c_ratio < 0.5:
        return False, "Oxygen-to-carbon ratio too low for lipopolysaccharide"
    
    # Check for presence of phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])[O-]")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    
    # LPS often contain phosphate groups
    if not has_phosphate:
        return False, "No phosphate groups found"
    
    # Check for presence of heptose and octulosonic acid units
    heptose_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@@H]([C@@H](CO)O1)O)O)O)O)O")
    octulosonic_acid_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@@H]([C@H](C(=O)O)O1)O)O)O)O")
    has_heptose = mol.HasSubstructMatch(heptose_pattern)
    has_octulosonic_acid = mol.HasSubstructMatch(octulosonic_acid_pattern)
    
    if not (has_heptose and has_octulosonic_acid):
        return False, "No heptose or octulosonic acid units found"
    
    # Check for presence of oligosaccharide chains
    oligosaccharide_pattern = Chem.MolFromSmarts("[OX2][CX4][OX2][CX4][OX2][CX4][OX2][CX4][OX2]")
    has_oligosaccharide = mol.HasSubstructMatch(oligosaccharide_pattern)
    
    if not has_oligosaccharide:
        return False, "No oligosaccharide chains found"
    
    # Check for presence of lipid chains
    lipid_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    has_lipid = mol.HasSubstructMatch(lipid_pattern)
    
    if not has_lipid:
        return False, "No lipid chains found"
    
    # Calculate topological fingerprint similarity to known LPS structures
    fps = [FingerprintMols.FingerprintMol(mol)]
    known_lps_mols = [Chem.MolFromSmiles(smi) for smi in known_lps_smiles]
    known_lps_fps = [FingerprintMols.FingerprintMol(m) for m in known_lps_mols]
    
    similarity_scores = [DataStructs.DiceSimilarity(fps[0], fp) for fp in known_lps_fps]
    avg_similarity = sum(similarity_scores) / len(similarity_scores)
    
    # LPS structures should have high similarity to known examples
    if avg_similarity < 0.6:
        return False, "Low structural similarity to known lipopolysaccharides"
    
    return True, "Molecule exhibits structural features characteristic of lipopolysaccharides"

# List of SMILES strings for known lipopolysaccharide structures
known_lps_smiles = [
    "O=C(OC(CCCC(=O)O[C@@H]1[C@@H](O)[C@H](O[C@@H]([C@H]1O)CO)O[C@H]2O[C@@H]([C@@H](O)[C@@H]([C@H]2O)O)CO)CCCCC)C(=O)O",
    "O=C(O[C@@H]1[C@@H](O[C@H](CO)[C@H]([C@@H]1OC(=O)CCCCCCCCCCCCCCC)OC(=O)C)OC[C@@H](O)[C@@H](O)CO)C(C)C",
    "CCCCCCCCCC(=O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O[C@@H](C[C@H](OC(C)=O)C(\C)=C\CO[C@H]([C@H](O)CO)[C@H](O)[C@H](O)CO)[C@](C)(O)CCC=C(C)C",
    "O=C(O[C@@H]1[C@@H](O[C@H](COC(=O)C)[C@H]([C@@H]1OC(=O)CCCCCCCCCCCCCCC)OC(=O)C)OC[C@@H](O)[C@@H](O)CO)CCCCC",
    "O=C(O[C@@H]1[C@@H](O[C@H](CO)[C@H]([C@@H]1OC(=O)CCCCCCCCCCCCC)OC(=O)C)OC[C@@H](O)[C@@H](O)CO)C(C)C"
]