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
    A sphingomyelin d18:1 has a sphingosine backbone (18 carbons, one double bond at position 4)
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

    # Check for phosphocholine group: O-P(=O)(O-)OCC[N+](C)(C)C
    phosphocholine_pattern = Chem.MolFromSmarts("[O][P](=[O])([O-])OCC[N+]($(C)(C)C)")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Missing phosphocholine group"

    # Check for amide-linked fatty acid: N-C(=O)
    amide_pattern = Chem.MolFromSmarts("[NH][C](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found"

    # Check sphingosine backbone: two hydroxyls and a trans double bond in the backbone
    # The core structure: [C@H](O)[C@@H](/C=C/CCCCCCCCCCCCC)N
    sphingosine_pattern = Chem.MolFromSmarts("[C@H]([OH])[C@@H](/C=C/CCCCCCCCCCCCC)N")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        # Try without stereochemistry
        sphingosine_pattern_fallback = Chem.MolFromSmarts("[CH]([OH])[CH](/C=C/CCCCCCCCCCCCC)N")
        if not mol.HasSubstructMatch(sphingosine_pattern_fallback):
            return False, "Sphingosine backbone not found"

    # Verify the amide is connected to the sphingosine's nitrogen
    # Get the nitrogen in the sphingosine pattern
    sphingosine_n = None
    for match in mol.GetSubstructMatches(sphingosine_pattern):
        sphingosine_n = match[2]  # Assuming the N is the third atom in the SMARTS
        break
    if not sphingosine_n:
        for match in mol.GetSubstructMatches(sphingosine_pattern_fallback):
            sphingosine_n = match[2]
            break

    if sphingosine_n is None:
        return False, "Sphingosine nitrogen not found"

    # Check if the sphingosine N is connected to the amide's C=O
    amide_n = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'N' and any(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in atom.GetBonds())]
    valid_amide = False
    for n_idx in amide_n:
        for bond in mol.GetAtomWithIdx(n_idx).GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                neighbor = bond.GetOtherAtomIdx(n_idx)
                if mol.GetAtomWithIdx(neighbor).GetSymbol() == 'C' and any(b.GetBondType() == Chem.rdchem.BondType.DOUBLE for b in mol.GetAtomWithIdx(neighbor).GetBonds()):
                    valid_amide = True
                    if n_idx == sphingosine_n:
                        break
    if not valid_amide:
        return False, "Amide not connected to sphingosine backbone"

    return True, "Sphingosine backbone with phosphocholine and amide-linked fatty acid"