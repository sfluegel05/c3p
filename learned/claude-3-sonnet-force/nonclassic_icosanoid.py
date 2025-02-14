"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: CHEBI:36738 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.DataStructs import BulkTanimotoSimilarity
import os

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string and fingerprint similarity.
    A nonclassic icosanoid is a biologically active signalling molecule made by oxygenation of C20 fatty acids
    other than the classic icosanoids (the leukotrienes and the prostanoids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check for C20 or fewer carbons
    if c_count > 20:
        return False, f"Molecule has more than 20 carbon atoms"
    
    # Look for carboxyl group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group (-COOH) found"
    
    # Look for at least 3 oxygens
    if o_count < 3:
        return False, "Fewer than 3 oxygen atoms found"
    
    # Check for absence of leukotriene substructure
    leukotriene_pattern = Chem.MolFromSmarts("CC(=O)C=CC=CC=CC=CC=CC=CC(=O)O")
    if mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Leukotriene substructure found, not a nonclassic icosanoid"
    
    # Check for absence of prostanoid substructure
    prostanoid_pattern = Chem.MolFromSmarts("CC(=O)C=CC=CC=CC=CC=CC=C")
    if mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Prostanoid substructure found, not a nonclassic icosanoid"

    # Load positive examples
    pos_examples = []
    pos_dir = os.path.join(os.path.dirname(__file__), "positive_examples")
    for filename in os.listdir(pos_dir):
        with open(os.path.join(pos_dir, filename), "r") as f:
            for line in f:
                pos_examples.append(Chem.MolFromSmiles(line.strip()))

    # Calculate fingerprint similarity to positive examples
    mol_fp = FingerprintMols.FingerprintMol(mol, FingerprintMols.PatternFingerprint())
    pos_fps = [FingerprintMols.FingerprintMol(ex, FingerprintMols.PatternFingerprint()) for ex in pos_examples]
    similarities = BulkTanimotoSimilarity(mol_fp, pos_fps)

    if max(similarities) >= 0.8:
        return True, f"Fingerprint similarity of {max(similarities):.2f} to known nonclassic icosanoid"
    else:
        return False, f"Fingerprint similarity of {max(similarities):.2f} to known nonclassic icosanoids is too low"