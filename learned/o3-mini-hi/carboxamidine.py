"""
Classifies: CHEBI:35359 carboxamidine
"""
#!/usr/bin/env python
"""
Classifies: Carboxamidine-containing compounds.
Definition: Compounds having the structure RC(=NR)NR2,
which is used in nomenclature to denote the -C(=NH)NH2 (or substituted analogues) group.
Examples include formamidine, acetamidine, benzamidine, guanidine, etc.
"""

from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine functional group (-C(=NH)NH2 or substituted analogues).

    The pattern looks for a trigonal carbon (CX3) that is double-bonded to an imino nitrogen (NX2)
    and single-bonded to an amine nitrogen (NX3), reflecting the RC(=NR)NR2 grouping.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a carboxamidine group, False otherwise.
        str: Explanation/reason of the classification.
    """
    # Parse the SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define the SMARTS pattern to find a carboxamidine group:
    # The pattern "[CX3](=[NX2])[NX3]" matches a trigonal carbon (C with 3 connections)
    # which is double bonded to an NX2 nitrogen and single bonded to an NX3 nitrogen.
    # This is intended to capture groups of the form RC(=NH)NH2 or its substituted derivatives.
    carboxamidine_pattern = Chem.MolFromSmarts("[CX3](=[NX2])[NX3]")
    if carboxamidine_pattern is None:
        return False, "Error defining the carboxamidine SMARTS pattern."

    # Use substructure matching.
    matches = mol.GetSubstructMatches(carboxamidine_pattern)
    if len(matches) == 0:
        return False, "Carboxamidine moiety (-C(=NH)NH2 or substituted equivalent) not found."
    
    # If at least one match is found, we assume the functional group is present.
    return True, f"Found carboxamidine group in {len(matches)} location(s)."

# For testing purposes:
if __name__ == '__main__':
    # List of example SMILES from the prompt (some examples):
    test_smiles = [
        "[H]C(N)=N",  # formamidine
        "CN(Cc1ccc(Cl)nc1)C(\\C)=N/C#N",  # (Z)-acetamiprid
        "CC(=N)NCC1=CC(CN)=CC=C1",  # N-[3-(aminomethyl)benzyl]acetamidine
        "NC(=N)c1ccccc1",  # benzamidine
        "NC(N)=N",  # guanidine
        "N=C(N)C=C"  # acrylamidine
    ]
    
    for smi in test_smiles:
        result, reason = is_carboxamidine(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")