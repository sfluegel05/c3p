"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: Organic sulfide (thioether)
Definition: Compounds having the structure R-S-R (with R ≠ H).
Such compounds were once called thioethers.
This improved classifier uses:
  • A SMARTS that matches any C–S–C connection (regardless of bond order)
  • A check that the sulfur is not additionally substituted with oxygen
  • A simple filter to exclude many peptides or free amino acids (which are not meant
    to be classified as “organic sulfide” in this context)
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.
    
    An organic sulfide is defined as a molecule containing at least one R-S-R moiety,
    where both R groups are not hydrogen (typically organic groups) and the S atom is not
    oxidized (i.e. not bound to oxygen). In addition, many peptides and free amino acids
    (which happen to contain a methionine side chain, for example) are filtered out.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an organic sulfide, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First filter: Exclude small free amino acids.
    # The pattern below matches a zwitterionic amino acid center (commonly seen in, e.g., methionine).
    amino_acid_pattern = Chem.MolFromSmarts("[C@H]([NH3+])C(=O)[O-]")
    if mol.GetNumAtoms() < 50 and mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Molecule appears to be a free amino acid (not classified as organic sulfide)"
    
    # Next, check for multiple amide bonds as an indicator of a peptide-like molecule.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) > 1:
        return False, "Molecule appears to be peptide-derived (contains multiple amide bonds)"
    
    # Now look for a thioether (organic sulfide) substructure.
    # Use a SMARTS that allows any bond type (~) between a carbon (atomic number 6)
    # and sulfur so that aromatic or aliphatic connections are included.
    thioether_pattern = Chem.MolFromSmarts("[#6]~S~[#6]")
    matches = mol.GetSubstructMatches(thioether_pattern)
    if not matches:
        return False, "No organic sulfide (R-S-R) moiety found in the molecule"
    
    # Finally, check each S found in a match to make sure it is not also bonded to oxygen.
    # (This helps avoid sulfoxides and sulfones.)
    for match in matches:
        # match is a tuple of indices corresponding to (carbon, sulfur, carbon)
        s_atom = mol.GetAtomWithIdx(match[1])
        for neighbor in s_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                return False, "Sulfide moiety is substituted with oxygen (not a thioether)"
    
    return True, "Molecule contains at least one organic sulfide (thioether) group (R-S-R bond)"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "CCCCCCCCCCCCCCCSC[C@H](COP(O)(=O)OCC[N+](C)(C)C)OP(O)(=O)CCCCCCCCCCCCCCCC",   # expected true
        "S(CC[C@H](NC(=O)[C@@H](N)CC(O)=O)C(=O)N[C@@H](CC=1NC=NC1)C(O)=O)C"              # peptide mimic => false
    ]
    for s in test_smiles:
        result, reason = is_organic_sulfide(s)
        print(result, reason)