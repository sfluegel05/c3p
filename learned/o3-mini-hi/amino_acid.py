"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies chemical entities of the class amino acid.
Definition used: A molecule is considered an amino acid if it contains at least one free carboxylic acid group
(and its deprotonated form) and at least one free (neutral) amino group (i.e. not acylated), and the two functional groups
occur “close” together (i.e. in the same amino acid residue). We measure closeness by requiring that at least one pair 
of acid (carboxyl carbon) and amino (nitrogen) atoms have a shortest bond path of 4 or fewer bonds.
This rule largely excludes di- or oligopeptides where the free acid and amino groups lie on different residues.
Examples that pass include (S)-gabaculine, robenacoxib, L-thyroxine and many other true amino acids.
Note: Some edge cases will remain ambiguous.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    Our algorithm requires:
      1. The molecule can be parsed.
      2. It contains at least one carboxylic acid group (either as a neutral acid [C(=O)[OH]] or an anion [C(=O)[O-]]).
      3. It contains at least one free (non-amidated, neutral) amino group.
      4. The free acid and amino groups occur “close together” (i.e. there is at least one pair
         for which the shortest bond path length is 4 bonds or fewer). This heuristic helps
         to avoid classifying peptides (where the free groups are on different amino acid units) as amino acids.
         
    Args:
         smiles (str): SMILES string of the molecule.
         
    Returns:
         bool: True if the molecule is classified as an amino acid, False otherwise.
         str: Explanation for the decision.
    """
    # Try to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define SMARTS for free (non-amidated) carboxylic acid groups.
    # Neutral acid: C(=O)[OH]
    acid_neutral_smarts = "[CX3](=O)[OX2H]"
    # Anionic acid: C(=O)[O-]
    acid_anion_smarts   = "[CX3](=O)[O-]"
    acid_neutral = Chem.MolFromSmarts(acid_neutral_smarts)
    acid_anion = Chem.MolFromSmarts(acid_anion_smarts)
    
    acid_matches = []
    for match in mol.GetSubstructMatches(acid_neutral):
        # In our pattern the acid carbon is the first atom in the match
        acid_matches.append(match[0])
    for match in mol.GetSubstructMatches(acid_anion):
        acid_matches.append(match[0])
        
    if not acid_matches:
        return False, "No carboxylic acid group found"
        
    # Define SMARTS for free amino groups.
    # We require the nitrogen to be sp3 (or not aromatic) and NOT directly attached to a carbonyl 
    # (i.e. not an amide) and also be neutral.
    # This pattern will match primary or secondary amines.
    amino_smarts = "[#7;X3;H1,H2;!$([#7]C=O);!$([#7+])]"
    amino = Chem.MolFromSmarts(amino_smarts)
    amino_matches = []
    for match in mol.GetSubstructMatches(amino):
        # In our pattern the nitrogen (free amine) is the first (and only) atom in the match
        amino_matches.append(match[0])
    if not amino_matches:
        return False, "No free (non-amidated) amino group found"
        
    # Now try to ensure that one acid and one amino group occur on the same “residue”
    # by calculating the shortest path (number of bonds) between any acid carbon (as identified above)
    # and any free amino nitrogen. In a typical amino acid these will be very close (2–4 bonds).
    # In peptides the free groups belong to different residues so the shortest path is considerably longer.
    min_dist = None
    for a_idx in acid_matches:
        for n_idx in amino_matches:
            path = rdmolops.GetShortestPath(mol, a_idx, n_idx)
            # path is a tuple of atom indices; length in bonds = len(path)-1
            if path:
                dist = len(path) - 1
                if min_dist is None or dist < min_dist:
                    min_dist = dist

    # If we were unable to get any reasonable path (should not happen for connected molecules)
    if min_dist is None:
        return False, "Could not determine connectivity between functional groups"
        
    # Use a threshold on bond distance.
    # (For a simple amino acid unit the free acid and amino group are typically separated by 2–4 bonds.
    # Peptides where the two free groups belong to different residues tend to have a longer distance.)
    cutoff = 4
    if min_dist <= cutoff:
        return True, f"Found free carboxylic acid and free amino group with a bond distance of {min_dist} (within cutoff)."
    else:
        return False, f"Functional groups found but they are separated by {min_dist} bonds, suggesting they belong to separate fragments (e.g. a peptide) rather than a single amino acid."

# Example test cases (can be removed in production code)
if __name__ == "__main__":
    test_cases = {
        "(S)-gabaculine": "C1=C(C[C@@H](C=C1)N)C(O)=O",
        "robenacoxib": "CCc1ccc(Nc2c(F)c(F)cc(F)c2F)c(CC(O)=O)c1",
        "His-Gly-Ser (dipeptide false positive)": "O=C(NCC(=O)N[C@@H](CO)C(O)=O)[C@@H](N)CC=1NC=NC1",
        "N-(methoxyacetyl)-4-hydroxyproline (missed amino group)": "COCC(=O)N1CC(O)CC1C(O)=O"
    }
    for name, sm in test_cases.items():
        result, reason = is_amino_acid(sm)
        print(f"{name}: {result} ({reason})")