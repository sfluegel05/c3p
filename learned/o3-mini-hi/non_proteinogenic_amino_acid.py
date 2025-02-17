"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
#!/usr/bin/env python
"""
Classifies: Non-proteinogenic amino acid 
Definition: Any amino-acid that is not naturally encoded in the genetic code of any organism.

The function is_non_proteinogenic_amino_acid takes a SMILES string as input and returns:
    (bool, str) -> (True, reason) if the molecule appears to be an amino acid that is NOT one of the canonical ones,
                  (False, reason) otherwise.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# Pre-compute InChIKeys for the 20 standard proteinogenic amino acids.
# (Glycine is included; note that in glycine there is no chirality.)
_PROTEINOGENIC_SMILES = [
    "NCC(=O)O",                         # Glycine
    "CC(N)C(=O)O",                       # Alanine
    "CC(C)[C@H](N)C(=O)O",                # Valine
    "CCC[C@H](N)C(=O)O",                  # Leucine
    "CC[C@H](C)[C@H](N)C(=O)O",           # Isoleucine
    "C1CC[C@@H](NC1)C(=O)O",              # Proline (L-Proline)
    "N[C@@H](Cc1ccccc1)C(=O)O",           # Phenylalanine
    "N[C@@H](Cc1ccc(O)cc1)C(=O)O",         # Tyrosine
    "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",    # Tryptophan
    "N[C@@H](CO)C(=O)O",                  # Serine
    "N[C@@H]([C@@H](O)C)C(=O)O",          # Threonine
    "N[C@@H](CS)C(=O)O",                  # Cysteine
    "N[C@@H](CCSC)C(=O)O",                # Methionine
    "N[C@@H](CC(=O)O)C(=O)O",             # Aspartic acid
    "N[C@@H](CCC(=O)O)C(=O)O",            # Glutamic acid
    "N[C@@H](CC(=O)N)C(=O)O",             # Asparagine
    "N[C@@H](CCC(=O)N)C(=O)O",            # Glutamine
    "N[C@@H](CCCCN)C(=O)O",               # Lysine
    "N[C@@H](CCCNC(=N)N)C(=O)O",          # Arginine
    "N[C@@H](Cc1cnc[nH]1)C(=O)O",          # Histidine
]

# Build a set of canonical InChIKeys for the proteinogenic amino acids.
_PROTEINOGENIC_INCHIKEYS = set()
for smi in _PROTEINOGENIC_SMILES:
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        # Use MolToInchiKey to compute a canonical identifier.
        ik = Chem.MolToInchiKey(mol)
        _PROTEINOGENIC_INCHIKEYS.add(ik)

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Classifies whether the given SMILES string represents a non-proteinogenic amino acid.
    
    The function first checks if the molecule contains an amino acid motif – that is, an alpha
    carbon bearing an amino group and a carboxylic acid group. Two SMARTS patterns are used:
        - One for amino acids with a chiral (substituted) alpha carbon.
        - One for glycine (which is achiral).
    
    Then, it computes the molecule's InChIKey and compares it against a set of canonical proteinogenic 
    amino acids. If a match is found the molecule is considered proteinogenic.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where bool is True if the molecule is a non-proteinogenic amino acid,
                     False otherwise; the reason is provided as a str.
    """
    
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the typical alpha-carbon motif in amino acids.
    # Pattern 1 matches amino acids with a chiral (substituted) alpha-carbon:
    #   An sp3 carbon with one hydrogen, attached to an amino group and a carboxyl group.
    patt_chiral = Chem.MolFromSmarts("[C;H1]([NX3])[C](=O)[O;H1,O-]")
    # Pattern 2 to match glycine (which has an unsubstituted alpha carbon with two hydrogens):
    patt_glycine = Chem.MolFromSmarts("NCC(=O)[O;H1,O-]")
    
    if not (mol.HasSubstructMatch(patt_chiral) or mol.HasSubstructMatch(patt_glycine)):
        return False, "Does not appear to contain the typical amino acid alpha-carbon motif"
    
    # Compute the InChIKey of the molecule.
    try:
        inchi_key = Chem.MolToInchiKey(mol)
    except Exception as e:
        return False, f"Failed to generate InChIKey: {e}"
    
    # Check if the input exactly matches one of the canonical proteinogenic amino acids.
    if inchi_key in _PROTEINOGENIC_INCHIKEYS:
        return False, "Matches a canonical proteinogenic amino acid"
    
    # If it passed both tests, we consider it an amino acid that is not encoded in the genetic code.
    return True, "Contains the amino acid motif and does not match any canonical proteinogenic amino acid"

# Example usage (uncomment for testing):
# test_smiles_list = [
#     "[C@H](N)(C(=O)O)",  # Glycine (should be proteinogenic)
#     "NCCCC[C@H](N)CC(O)=O",  # (3S)-3,7-diaminoheptanoic acid (non-proteinogenic)
# ]
# for smi in test_smiles_list:
#     result, reason = is_non_proteinogenic_amino_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")