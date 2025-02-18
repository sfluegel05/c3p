"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: proteinogenic amino acid
Defined as one of the 23 α‐amino acids that are incorporated into proteins.
Examples include L-glutamic acid, L-histidine, glycine (non-chiral), proline (cyclic) etc.
"""

from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    The heuristic is based on the observation that every proteinogenic amino acid has
    an α‐carbon attached to both a carboxyl group and an amino group. For glycine no
    chiral tag is present, and for proline the amino group is contained in a ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a proteinogenic amino acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns that capture common motifs for proteinogenic amino acids.
    # Pattern 1 and 2: A chiral α-carbon with an amino group and a carboxyl group.
    patterns = [
        "[C@H](N)C(=O)[O]",   # chiral center with amino and carboxyl (one stereochemistry)
        "[C@@H](N)C(=O)[O]",  # mirror image stereochemistry
        # Pattern for glycine (no chiral tag)
        "NC(C(=O)[O])",
        # Patterns capturing proline-type amino acids (the amino group is in a ring)
        "[C@H]1CCCN1C(=O)[O]",
        "[C@@H]1CCCN1C(=O)[O]"
    ]
    
    # Compile SMARTS patterns.
    smarts_mols = []
    for smarts in patterns:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            smarts_mols.append((smarts, patt))
    
    # Check if any SMARTS pattern is found as a substructure.
    for smarts, patt in smarts_mols:
        if mol.HasSubstructMatch(patt):
            return True, f"Matches the pattern: {smarts}"
    
    # If no patterns match, the molecule does not appear to be a proteinogenic amino acid.
    return False, "No appropriate α-amino acid motif (amino and carboxyl groups on the same carbon) found"

# For quick ad hoc testing, one might run:
if __name__ == "__main__":
    # Examples provided (picked two examples: L-glutamic acid and proline)
    test_smiles = {
        "L-glutamic acid": "N[C@@H](CCC(O)=O)C(O)=O",
        "L-proline": "OC(=O)[C@@H]1CCCN1",
        "Invalid": "NotAValidSMILES"
    }
    for name, smi in test_smiles.items():
        result, reason = is_proteinogenic_amino_acid(smi)
        print(f"{name}: {result} ({reason})")