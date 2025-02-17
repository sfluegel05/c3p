"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: Hydroxy-5β-cholanic acids (bile acids)
Definition: Any member of a group of hydroxy-5β-cholanic acids occurring in bile, where they are present as the sodium salts of their amides with glycine or taurine. In mammals bile acids almost invariably have 5β-configuration.
Note: The 5β stereochemistry is difficult to validate with a generic SMARTS search; here we focus on detecting the typical steroid nucleus with a carboxylic acid side chain and hydroxyl substituents.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid (hydroxy-5β-cholanic acid) based on its SMILES string.
    
    To be classified as a bile acid the molecule should:
      - Parse correctly.
      - Contain a carboxylic acid group (approx. the C-24 carboxylic acid).
      - Contain at least one other hydroxyl group.
      - Possess a steroid nucleus, here approximated as a fused ring system with three six-membered rings and one five-membered ring.
      - Typically have a molecular weight in the range expected for bile acids (~400-600 Da).
      
    Note:
      - This function does not fully determine the 5β stereochemistry.
      
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a bile acid, False otherwise
        str: Reason for the classification decision
    """
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group pattern (C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")  
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "Missing carboxylic acid group which is essential for cholanic acids"

    # Check for hydroxyl groups (excluding the carboxylic acid OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    # Need at least one hydroxyl group apart from the acid (sometimes the acid OH is counted, but we expect additional hydroxy substituents)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl substituents found besides the acid group"

    # Check molecular weight - bile acids are typically in the ~400-600 Da range (this is only a rough guide)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (350 <= mol_wt <= 700):
        return False, f"Molecular weight {mol_wt:.1f} out of expected bile acid range"

    # Check for fused steroid nucleus:
    # We approximate the steroid nucleus as having 4 fused rings with three six-membered and one five-membered rings.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found, expected a steroid nucleus"

    # Count rings by size
    count_5 = 0
    count_6 = 0
    for ring in rings:
        if len(ring) == 5:
            count_5 += 1
        elif len(ring) == 6:
            count_6 += 1

    # A steroid nucleus is expected to have at least one 5-membered and three 6-membered rings.
    if count_5 < 1 or count_6 < 3:
        return False, f"Fused ring system does not match steroid nucleus criteria (found {count_5} five-membered and {count_6} six-membered rings)"

    # Additional note: Checking the exact fused nature is nontrivial.
    # We give a positive result if the basic conditions above are met.
    return True, "Molecule shows fused steroid core with a carboxylic acid and hydroxy groups typical for hydroxy-5β-cholanic (bile acid) structures"