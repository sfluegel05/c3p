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
    
    # More general naphthoquinone patterns that capture both ortho and para arrangements
    # Using aromatic and non-aromatic carbons to be more flexible
    naphthoquinone_patterns = [
        # Basic naphthoquinone core with flexible bond types
        "[#6]1[#6]~[#6]~[#6]2[#6](=[O])[#6]~[#6](=[O])[#6]~[#6]12",
        # Alternative pattern with different quinone positions
        "[#6]1[#6]~[#6]~[#6]2[#6](=[O])[#6]~[#6]~[#6](=[O])[#6]12",
        # Pattern for more complex fused systems
        "[#6]2~[#6]~[#6]~[#6]1[#6](~[#6]~[#6](~[#6]2)~[#6](=[O])[#6]1)=[O]"
    ]
    
    has_naphthoquinone = False
    core_atoms = set()
    
    for pattern in naphthoquinone_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            has_naphthoquinone = True
            core_atoms.update(mol.GetSubstructMatch(patt))
    
    if not has_naphthoquinone:
        return False, "No naphthoquinone core found"
    
    # Look for both hydroxyl groups and deprotonated oxygens attached to the core
    hydroxyl_patterns = [
        "[OX2H1]",  # regular hydroxyl
        "[O-]",     # deprotonated hydroxyl
        "[OX2H0]-[#6]"  # other oxygen attachments that might be hydroxyls
    ]
    
    for pattern in hydroxyl_patterns:
        hydroxyl_pattern = Chem.MolFromSmarts(pattern)
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        
        for match in hydroxyl_matches:
            oh_oxygen = mol.GetAtomWithIdx(match[0])
            for neighbor in oh_oxygen.GetNeighbors():
                if neighbor.GetIdx() in core_atoms:
                    return True, "Contains naphthoquinone core with at least one hydroxyl group attached"
    
    return False, "No hydroxyl groups attached to naphthoquinone core"

def test_examples():
    """Test function with some example molecules"""
    examples = [
        ("Oc1cc(O)c2C(=O)C=C(O)C(=O)c2c1", "flaviolin"),
        ("OC1=C([C@H]2CC[C@@H](CC2)c2ccc(Cl)cc2)C(=O)c2ccccc2C1=O", "atovaquone"),
        ("CC(C)=CC[C@H](O)C1=CC(=O)c2c(O)ccc(O)c2C1=O", "Alkannin"),
        ("CC1=CC(=O)C=CC1=O", "not a hydroxynaphthoquinone (methylquinone)"),
        ("[C@H](C)([C@@H]([C@@H]([C@H](\C=C\O[C@]1(OC=2C(C1=O)=C3C(C(C(C(=C3[O-])/C=N/N4CCN(CC4)C)=O)=O)=C(C2C)[O-])C)OC)C)OC(=O)C)[C@H](O)[C@@H]([C@@H](O)[C@@H](C)/C=C/C=C(/C)\C(N)=O)C", "rifampicin derivative")
    ]
    
    for smi, name in examples:
        result, reason = is_hydroxynaphthoquinone(smi)
        print(f"Name: {name}")
        print(f"SMILES: {smi}")
        print(f"Is hydroxynaphthoquinone: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()