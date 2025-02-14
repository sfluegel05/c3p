"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthaoquinone moiety is substituted by at least one hydroxy group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for naphthoquinone cores
    # 1,4-naphthoquinone pattern
    pattern_14 = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]1):[#6]:[#6]C(=O)[#6]2=O")
    # 1,2-naphthoquinone pattern
    pattern_12 = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]1)C(=O)C(=O)[#6]2")
    
    if pattern_14 is None or pattern_12 is None:
        return None, "Error in SMARTS pattern"

    has_naphthoquinone = False
    core_atoms = set()
    
    # Check for 1,4-naphthoquinone
    if mol.HasSubstructMatch(pattern_14):
        has_naphthoquinone = True
        matches = mol.GetSubstructMatches(pattern_14)
        for match in matches:
            core_atoms.update(match)
    
    # Check for 1,2-naphthoquinone
    if mol.HasSubstructMatch(pattern_12):
        has_naphthoquinone = True
        matches = mol.GetSubstructMatches(pattern_12)
        for match in matches:
            core_atoms.update(match)

    if not has_naphthoquinone:
        return False, "No naphthoquinone core found"

    # Look for hydroxy groups
    # Count oxygens attached to the core
    hydroxy_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # oxygen atom
            # Get neighboring carbons
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetIdx() in core_atoms:
                    hydroxy_count += 1
                    break

    if hydroxy_count == 0:
        return False, "No hydroxy groups attached to naphthoquinone core"
    
    return True, f"Contains naphthoquinone core with {hydroxy_count} hydroxy group(s)"

def test_examples():
    """Test function with some example molecules"""
    examples = [
        ("Oc1cc(O)c2C(=O)C=C(O)C(=O)c2c1", "flaviolin"),
        ("OC1=C([C@H]2CC[C@@H](CC2)c2ccc(Cl)cc2)C(=O)c2ccccc2C1=O", "atovaquone"),
        ("CC(C)=CC[C@H](O)C1=CC(=O)c2c(O)ccc(O)c2C1=O", "Alkannin"),
        ("CC1=CC(=O)C=CC1=O", "methylquinone (negative)"),
        ("Oc1cccc2C(=O)C=CC(=O)c12", "juglone")
    ]
    
    for smi, name in examples:
        result, reason = is_hydroxynaphthoquinone(smi)
        print(f"Name: {name}")
        print(f"SMILES: {smi}")
        print(f"Is hydroxynaphthoquinone: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()