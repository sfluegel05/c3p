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
        
    # 2. Check for core porphyrin system with Mg (specific pattern)
    # The pattern ensures alternating double and single bonds, and the Mg coordination to 4 nitrogens
    porphyrin_core_pattern = Chem.MolFromSmarts("[n1X2][cX3][cX3][n2X2][c:1][c:2][n3X2][c:3][c:4][n4X2][Mg]1234")
    if not mol.HasSubstructMatch(porphyrin_core_pattern):
        return False, "No porphyrin core with magnesium found"
    
    # 3. Check for a fifth ring attached to a specific porphyrin carbon.
    #   The 5th ring is specifically attached to one of the carbon atoms of the porphyrin ring.
    #   The [c:1] maps to the core porphyrin's carbon ring system.
    fifth_ring_pattern = Chem.MolFromSmarts("[c:1]1[c][c][c][c]1")
    if not mol.HasSubstructMatch(fifth_ring_pattern):
        return False, "No fifth ring attached to the porphyrin core found"

    # 4. Check for at least one long phytol chain esterified to carboxyl
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if not ester_matches:
       return False, f"No ester group found"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    
    # Iterate through ester groups and check attached chains
    has_phytol_chain = False
    for match in ester_matches:
        ester_oxygen_atom = mol.GetAtomWithIdx(match[0]) # get the oxygen atom of the ester
        
        for neighbor in ester_oxygen_atom.GetNeighbors():
           if mol.HasSubstructMatch(fatty_acid_pattern, fromAtomIdx=neighbor.GetIdx()):
              has_phytol_chain = True
              break
        if has_phytol_chain:
            break    
    if not has_phytol_chain:
       return False, f"Missing phytol-like chain attached to the ester"

    return True, "Contains a magnesium porphyrin core, a fifth ring, and a phytol chain"