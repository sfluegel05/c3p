"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: L-alpha-amino acid
Definition: Any alpha-amino acid having L-configuration (that is, the α-carbon has 
the “S” CIP configuration) and a free (non-amidated, protonated) carboxylic acid group.

The function is_L_alpha_amino_acid determines if a molecule (given as SMILES) has exactly one
free amino acid backbone: a free (non‐amidated) NH2 group directly bonded to a chiral α‐carbon 
which is bonded to a carboxylic acid group (C(=O)[OH]). This function uses a pair of SMARTS patterns 
to detect the proper chirality (either [C@H] or [C@@H]) and then verifies that:
  1. The amino nitrogen is free (bound only to H’s and the α‑carbon).
  2. The corresponding carboxyl carbon is bonded to exactly two oxygen atoms (one double‐bonded,
     and one single‐bonded, ideally carrying a hydrogen).
  3. The α‑carbon’s CIP code is “S”.
If any of these checks fail, a meaningful error message is returned.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is a free L-alpha-amino acid.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a free L-alpha-amino acid, False otherwise.
        str: Reason explaining the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to inspect connectivity reliably.
    mol = Chem.AddHs(mol)
    
    # Compute stereochemistry and assign CIP codes.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define SMARTS patterns for the free amino acid backbone.
    # Looking for an amino nitrogen (not amidated) attached to a chiral α‐carbon bonded to a carbon side chain,
    # which in turn is bonded to a carboxyl carbon with an (assumed) free carboxylic acid.
    pattern1 = Chem.MolFromSmarts("[NX3;H1,H2;!$(N-C(=O))]-[C@H]([#6])-[C](=O)[O;H1]")
    pattern2 = Chem.MolFromSmarts("[NX3;H1,H2;!$(N-C(=O))]-[C@@H]([#6])-[C](=O)[O;H1]")
    
    matches1 = mol.GetSubstructMatches(pattern1, useChirality=True)
    matches2 = mol.GetSubstructMatches(pattern2, useChirality=True)
    
    total_matches = []
    for m in matches1:
        total_matches.append(("pattern1", m))
    for m in matches2:
        total_matches.append(("pattern2", m))
    
    # Remove duplicate matches (based on atom indices in the match)
    unique_matches = []
    seen = set()
    for label, match in total_matches:
        tup = tuple(match)
        if tup not in seen:
            seen.add(tup)
            unique_matches.append((label, match))
    
    if len(unique_matches) == 0:
        return False, "Alpha-amino acid backbone (free NH2 and COOH group) not found"
    if len(unique_matches) > 1:
        return False, f"Found {len(unique_matches)} amino acid backbone motifs; likely a peptide or more than one motif exists"
    
    # Our single match is our candidate backbone.
    label, match = unique_matches[0]
    # According to our SMARTS, atom indices are: 0: amino nitrogen, 1: α-carbon, 2: carboxyl carbon.
    N_idx, Ca_idx, Cc_idx = match[0], match[1], match[2]
    N_atom = mol.GetAtomWithIdx(N_idx)
    Ca_atom = mol.GetAtomWithIdx(Ca_idx)
    Cc_atom = mol.GetAtomWithIdx(Cc_idx)
    
    # Check that the amino nitrogen is free (only attached to the α-carbon as heavy neighbor).
    heavy_neigh_N = [nbr for nbr in N_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    if len(heavy_neigh_N) != 1 or heavy_neigh_N[0].GetIdx() != Ca_idx:
        return False, "Backbone nitrogen appears modified or is not free (unexpected heavy atom connectivity)"
    
    # Instead of enforcing that the carboxyl carbon has a strict degree of 3,
    # we now check that it is attached to exactly two oxygen atoms (aside from the α-carbon).
    oxy_neighbors = [nbr for nbr in Cc_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(oxy_neighbors) != 2:
        return False, "Carboxyl carbon does not have exactly two oxygen neighbors"
    
    # Verify that one oxygen is double-bonded (the carbonyl oxygen) and 
    # the other is single-bonded (the hydroxyl oxygen, preferably with at least one hydrogen).
    found_double = False
    found_single = False
    for o in oxy_neighbors:
        bond = mol.GetBondBetweenAtoms(Cc_atom.GetIdx(), o.GetIdx())
        if bond is None:
            continue
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            found_double = True
        elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            # Check if the single-bonded oxygen carries any hydrogens.
            if any(n.GetAtomicNum() == 1 for n in o.GetNeighbors()):
                found_single = True
    if not (found_double and found_single):
        return False, "The carboxylic acid group does not appear protonated/free (expected C(=O)[OH])"
    
    # Check the chirality of the α-carbon using the assigned CIP code.
    if not Ca_atom.HasProp("_CIPCode"):
        return False, "Alpha-carbon lacks a CIP code; cannot determine configuration"
    cip = Ca_atom.GetProp("_CIPCode")
    if cip != "S":
        return False, f"Alpha-amino acid backbone found but alpha-carbon CIP code is '{cip}', not 'S'"
    
    return True, "Found a free L-alpha-amino acid backbone (free NH2 and free COOH) with L (S) configuration"

# Example usage:
if __name__ == '__main__':
    # Testing with a few examples from the provided list
    test_smiles = [
        "N[C@@H](CC(=O)c1cc(O)ccc1N)C(O)=O",         # 5-hydroxy-L-kynurenine
        "N[C@@H](CC(=C)C(N)=O)C(O)=O",                # 4-methylene-L-glutamine
        "CC[C@H](N)C(O)=O",                          # L-alpha-aminobutyric acid
        "O=C(O)[C@@H](N)CCC1C=CC(N)C=C1",             # Amiclenomycin
        "N[C@@H](CCCCC(O)=O)C(O)=O",                  # L-2-aminopimelic acid
        "N[C@@H](CCCCNC(O)=O)C(O)=O",                 # N(6)-carboxy-L-lysine
        "CNC(=O)C[C@H](N)C(O)=O",                     # N(4)-methyl-L-asparagine
        "N[C@@H](COS(O)(=O)=O)C(O)=O",                # L-serine O-sulfate
        "N[C@@H](CCC(N)=O)C(O)=O",                    # L-glutamine
        "Cn1cncc1C[C@H](N)C(O)=O"                     # N(pros)-methyl-L-histidine
    ]
    for smi in test_smiles:
        res, reason = is_L_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*60}")