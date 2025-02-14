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

    # Define patterns for 1,2- and 1,4-naphthoquinone cores
    naphthoquinone_patterns = [
        # 1,4-naphthoquinone core
        "[#6]1:,[#6]-,:[#6]:,[#6]:,[#6]2:,[#6]:,[#6]:,[#6]:,[#6]:,[#6]1[#6](=O)-,:[#6]:,[#6]-,:[#6]2=O",
        # 1,2-naphthoquinone core
        "[#6]1:,[#6]-,:[#6]:,[#6]:,[#6]2:,[#6]:,[#6]:,[#6]:,[#6]:,[#6]1[#6](=O)-,:[#6]2=O"
    ]
    
    has_naphthoquinone = False
    matched_atoms = set()
    
    for pattern in naphthoquinone_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            has_naphthoquinone = True
            # Get all matches to find the core atoms
            matches = mol.GetSubstructMatches(patt)
            for match in matches:
                matched_atoms.update(match)

    if not has_naphthoquinone:
        return False, "No naphthoquinone core found"

    # Look for hydroxyl groups attached to the core
    # Include both -OH and -O- forms
    hydroxyl_patterns = [
        "[OX2H]-[#6]",  # neutral hydroxyl
        "[O-]-[#6]"     # deprotonated hydroxyl
    ]
    
    for pattern in hydroxyl_patterns:
        hydroxyl_patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(hydroxyl_patt):
            matches = mol.GetSubstructMatches(hydroxyl_patt)
            for match in matches:
                # Get the carbon atom the hydroxyl is attached to
                carbon_idx = match[1]  # second atom in the match is the carbon
                if carbon_idx in matched_atoms:
                    return True, "Contains naphthoquinone core with hydroxyl group attached"

    return False, "No hydroxyl groups attached to naphthoquinone core"

def test_examples():
    """Test function with some example molecules"""
    examples = [
        ("Oc1cc(O)c2C(=O)C=C(O)C(=O)c2c1", "flaviolin"),
        ("OC1=C([C@H]2CC[C@@H](CC2)c2ccc(Cl)cc2)C(=O)c2ccccc2C1=O", "atovaquone"),
        ("CC(C)=CC[C@H](O)C1=CC(=O)c2c(O)ccc(O)c2C1=O", "Alkannin"),
        ("CC1=CC(=O)C=CC1=O", "methylquinone (negative example)"),
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