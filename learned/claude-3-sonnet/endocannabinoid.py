"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: CHEBI:38608 endocannabinoid
"""
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
import os

# Load a database of known endocannabinoid SMILES strings
known_endocannabinoids = set()
with open("endocannabinoid_database.smi", "r") as f:
    for line in f:
        known_endocannabinoids.add(line.strip())

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Uses a database of known endocannabinoid structures and molecular similarity measures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if SMILES matches a known endocannabinoid
    if smiles in known_endocannabinoids:
        return True, "Matches a known endocannabinoid structure"

    # Compute molecular fingerprints
    fp_mol = MACCSkeys.GenMACCSKeys(mol)
    
    # Compute similarity to known endocannabinoids
    max_similarity = 0
    for endocannabinoid in known_endocannabinoids:
        endocannabinoid_mol = Chem.MolFromSmiles(endocannabinoid)
        fp_endocannabinoid = MACCSkeys.GenMACCSKeys(endocannabinoid_mol)
        similarity = DataStructs.FingerprintSimilarity(fp_mol, fp_endocannabinoid)
        max_similarity = max(max_similarity, similarity)

    # Classify based on similarity threshold
    similarity_threshold = 0.7
    if max_similarity >= similarity_threshold:
        return True, f"Structurally similar to known endocannabinoids (similarity: {max_similarity:.2f})"
    else:
        return False, f"Not structurally similar to known endocannabinoids (max similarity: {max_similarity:.2f})"