"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:17578 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    Must contain:
    1. Sphingosine backbone with amino alcohol and double bond
    2. Amide-linked fatty acid chain
    3. Minimum 12 carbons in sphingosine chain
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for amide group (N connected to carbonyl)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Find sphingosine backbone pattern: amino alcohol with double bond
    # Pattern: NH connected to CH(OH)CH2... with a double bond in the chain
    sphingosine_pattern = Chem.MolFromSmarts("[NH][C]([OH])([CH2])[CH2]/[CH]=[CH]/[CH2]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Missing sphingosine backbone features"

    # Verify fatty acid chain length (amide-connected chain)
    # Get amide nitrogen and follow R-group
    matches = mol.GetSubstructMatches(amide_pattern)
    for match in matches:
        nitrogen_idx = match[0]
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Get connected carbon chain (fatty acid)
        fatty_acid_chain = []
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != match[1]:
                chain_atom = neighbor
                break
        else:
            continue
        
        # Traverse the fatty acid chain
        chain_length = 0
        stack = [(chain_atom, 0)]
        visited = set()
        while stack:
            atom, depth = stack.pop()
            if atom.GetAtomicNum() != 6 or atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            chain_length = max(chain_length, depth + 1)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                    stack.append((nbr, depth + 1))
        
        if chain_length >= 8:  # Minimum 8 carbons in fatty acid
            return True, "Contains sphingosine backbone with N-linked fatty acyl chain"

    return False, "Fatty acid chain too short or missing required features"