"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are magnesium porphyrins with a fifth ring and usually a long phytol chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for magnesium
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found."
    if len(mg_atoms) != 1:
        return False, "Incorrect number of magnesium atoms."
        
    # 2. Check for core porphyrin system with Mg (relaxed pattern)
    # Look for a 4 nitrogen ring system coordinating to Mg and containing carbons.
    porphyrin_core_pattern = Chem.MolFromSmarts("[nX2,nX3][cX3][cX3][nX2,nX3][cX3][cX3][nX2,nX3][cX3][cX3][nX2,nX3][Mg]")
    if not mol.HasSubstructMatch(porphyrin_core_pattern):
        return False, "No porphyrin core with magnesium found"
    
    # 3. Check for a fifth ring attachment (relaxed pattern)
    # check for a ring atom fused to another ring.
    fifth_ring_pattern = Chem.MolFromSmarts("[R]1~[R]~[R]~[R]~[R]~1")
    if not mol.HasSubstructMatch(fifth_ring_pattern):
        return False, "No fifth ring attached to the porphyrin core found"
    
    # 4. Check for at least one long phytol chain esterified to carboxyl
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if not ester_matches:
       return False, f"No ester group found"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
       return False, f"Missing fatty acid chains"

    return True, "Contains a magnesium porphyrin core, a fifth ring, and a phytol chain"