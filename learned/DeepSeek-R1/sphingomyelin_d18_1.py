"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: CHEBI: sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    A sphingomyelin d18:1 has a sphingosine backbone (18 carbons, one trans double bond at position 4)
    with a phosphocholine group and an amide-linked fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin d18:1, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for phosphocholine group: O-P(=O)([O-])OCC[N+](C)(C)C
    phosphocholine_pattern = Chem.MolFromSmarts("[O]P(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Missing phosphocholine group"

    # Check sphingosine backbone: two hydroxyls and a trans double bond in the backbone
    # Core structure: [C@...](O)[C@@...](/C=C/CCCCCCCCCCCCC)N
    sphingosine_pattern = Chem.MolFromSmarts("[C@H]([OH])[C@@H](/C=C/CCCCCCCCCCCCC)N")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        # Try without stereochemistry but require /C=C/ (trans)
        sphingosine_pattern_fallback = Chem.MolFromSmarts("[CH]([OH])[CH](/C=C/CCCCCCCCCCCCC)N")
        if not mol.HasSubstructMatch(sphingosine_pattern_fallback):
            return False, "Sphingosine backbone not found"

    # Find the sphingosine nitrogen atom
    sphingosine_n = None
    for n_pos in mol.GetSubstructMatches(sphingosine_pattern):
        sphingosine_n = n_pos[-1]  # N is last atom in SMARTS pattern
        break
    if sphingosine_n is None:
        for n_pos in mol.GetSubstructMatches(sphingosine_pattern_fallback):
            sphingosine_n = n_pos[-1]
            break
    if sphingosine_n is None:
        return False, "Sphingosine nitrogen not found"

    # Verify the sphingosine N is part of an amide (N connected to C=O)
    n_atom = mol.GetAtomWithIdx(sphingosine_n)
    has_amide = False
    for bond in n_atom.GetBonds():
        neighbor = bond.GetOtherAtomIdx(sphingosine_n)
        neighbor_atom = mol.GetAtomWithIdx(neighbor)
        if neighbor_atom.GetSymbol() == 'C':
            # Check if this C is double-bonded to an O
            for b in neighbor_atom.GetBonds():
                if b.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    other_atom = b.GetOtherAtom(neighbor_atom)
                    if other_atom.GetSymbol() == 'O':
                        has_amide = True
                        break
            if has_amide:
                break
    if not has_amide:
        return False, "Sphingosine nitrogen not part of amide bond"

    return True, "Sphingosine backbone with phosphocholine and amide-linked fatty acid"