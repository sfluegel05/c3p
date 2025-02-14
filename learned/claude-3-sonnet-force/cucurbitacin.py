"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from the triterpene hydrocarbon cucurbitane,
    characterized by a rigid tetracyclic scaffold with oxygenated substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cucurbitane core
    cucurbitane_pattern = Chem.MolFromSmarts("[C@@]12[C@H]([C@@]3([C@H]([C@@H]4[C@@H](C(=C3)C)C[C@@H](C4)O)CC2)C)[C@H](C1)O")
    if not mol.HasSubstructMatch(cucurbitane_pattern):
        return False, "No cucurbitane core found"

    # Check for oxygenated substituents (hydroxy, keto, ester, ether, glycosides)
    oxy_patterns = [Chem.MolFromSmarts("[OX2]"),
                    Chem.MolFromSmarts("[OX1]=C"),
                    Chem.MolFromSmarts("[OX2][CX3](=O)"),
                    Chem.MolFromSmarts("[OX2][CX3]")]
    oxy_matches = sum(mol.HasSubstructMatch(pat) for pat in oxy_patterns)
    if oxy_matches < 2:
        return False, "Not enough oxygenated substituents found"

    # Check molecular weight - cucurbitacins typically 500-1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for cucurbitacins"

    # Check for tetracyclic scaffold
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    if num_rings < 4:
        return False, "Not a tetracyclic compound"

    # Count carbons, oxygens, and rings
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20 or c_count > 40:
        return False, "Carbon count outside typical range for cucurbitacins"
    if o_count < 3 or o_count > 10:
        return False, "Oxygen count outside typical range for cucurbitacins"

    return True, "Contains cucurbitane core with oxygenated substituents and tetracyclic scaffold"