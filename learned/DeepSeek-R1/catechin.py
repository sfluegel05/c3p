"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are hydroxyflavans with a flavan-3-ol skeleton and substituted derivatives.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Core flavan-3-ol pattern: benzopyran system with 3-OH
    # Matches O1C2=CC=CC=C2C[C@H](O)C1 (chromanol structure)
    flavanol_core = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]-2-[#6][C@H](O)[#6][O;H1]1>>2")
    if not mol.HasSubstructMatch(flavanol_core):
        return False, "Missing flavan-3-ol core (chromanol with 3-OH)"

    # Verify two aromatic rings (A and B rings) in the chromane system
    ring_info = mol.GetRingInfo()
    aromatic_rings = [ring for ring in ring_info.AtomRings() 
                     if len(ring) >= 6 and all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)]
    
    # Check if core connects two aromatic rings
    core_atoms = set(flavanol_core.GetSubstructMatch(mol))
    connected_aromatics = 0
    for ring in aromatic_rings:
        if any(a in core_atoms for a in ring):
            connected_aromatics += 1
    
    if connected_aromatics < 2:
        return False, f"Only {connected_aromatics} aromatic rings connected to core"

    # Allow derivatives by checking for at least 2 oxygen atoms (core OH + substituents)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:
        return False, f"Only {oxygen_count} oxygen atoms (need ≥2)"

    return True, "Contains flavan-3-ol core with two connected aromatic rings"