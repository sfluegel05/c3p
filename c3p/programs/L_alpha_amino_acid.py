"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: L-alpha-amino acid
Definition: A free L-alpha-amino acid is defined as a molecule having one (and only one)
free amino acid backbone. That is, it contains an unmodified (non-amidated) amino group
(NH2) attached to a chiral α‐carbon which is bonded to a protonated carboxylic acid 
group (–C(=O)[OH]) and the α‑carbon has the “S” (L) CIP configuration.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is a free L-alpha-amino acid.
    Checks:
      1. That the molecule parses and hydrogens are added.
      2. That exactly one free amino acid backbone (NH2–C*–COOH) is found.
         This is done with two SMARTS patterns (one for each chiral annotation).
      3. The backbone is “free” – that is, the amino nitrogen is only bound to the α‐carbon, 
         and the carboxyl carbon is only bound to the α‐carbon and exactly two oxygens,
         one of which is single‐bonded and bears a hydrogen, and the other double‐bonded.
      4. That the α‐carbon has a CIP code (assigned by RDKit) and that it is “S”.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a free L-alpha-amino acid, False otherwise.
        str: Reason explaining the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens so that we can inspect the connectivity and H counts.
    mol = Chem.AddHs(mol)
    
    # Compute stereochemistry and assign CIP codes.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define SMARTS patterns for the free amino acid backbone.
    # These patterns look for an amino nitrogen that is not amidated,
    # attached to a chiral α‐carbon (either @ or @@) bearing a carbon side chain ([#6])
    # and then attached to a carboxyl group of the form C(=O)[O;H1].
    pattern1 = Chem.MolFromSmarts("[NX3;H1,H2;!$(N-C(=O))]-[C@H]([#6])-[C](=O)[O;H1]")
    pattern2 = Chem.MolFromSmarts("[NX3;H1,H2;!$(N-C(=O))]-[C@@H]([#6])-[C](=O)[O;H1]")
    
    matches1 = mol.GetSubstructMatches(pattern1, useChirality=True)
    matches2 = mol.GetSubstructMatches(pattern2, useChirality=True)
    
    total_matches = []
    for m in matches1:
        total_matches.append(("pattern1", m))
    for m in matches2:
        total_matches.append(("pattern2", m))
    
    # Use unique matches based on the raw indices from the match, but preserve the ordering from match.
    unique_matches = []
    seen = set()
    for label, match in total_matches:
        tup = tuple(match)
        if tup not in seen:
            seen.add(tup)
            unique_matches.append((label, match))
    
    if len(unique_matches) == 0:
        return False, "Alpha-amino acid backbone (free NH2 and protonated COOH) not found"
    if len(unique_matches) > 1:
        return False, f"Found {len(unique_matches)} amino acid backbones; likely a peptide or multiple motifs exist"
    
    # Get the single match.
    label, match = unique_matches[0]
    # According to our SMARTS, the ordering is: 0: amino nitrogen, 1: α-carbon, 2: carboxyl carbon.
    N_idx, Ca_idx, Cc_idx = match[0], match[1], match[2]
    N_atom = mol.GetAtomWithIdx(N_idx)
    Ca_atom = mol.GetAtomWithIdx(Ca_idx)
    Cc_atom = mol.GetAtomWithIdx(Cc_idx)
    
    # --- Check that the backbone is free (unmodified) ---
    # For the amino nitrogen: in a free amino acid, it should be NH2 (only bound to hydrogens besides the α‑carbon).
    heavy_neigh_N = [nbr for nbr in N_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    if len(heavy_neigh_N) != 1 or heavy_neigh_N[0].GetIdx() != Ca_idx:
        return False, "Backbone nitrogen appears modified or not free (unexpected heavy atom connectivity)"
    
    # For the carboxyl carbon: It should be bound to exactly three atoms:
    # the α‑carbon and two oxygens.
    if Cc_atom.GetDegree() != 3:
        return False, "Carboxyl carbon connectivity is not as expected for a free COOH group"
    # Identify the two oxygens attached to Cc_atom.
    oxy_atoms = [nbr for nbr in Cc_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(oxy_atoms) != 2:
        return False, "Carboxyl carbon does not have exactly two oxygen neighbors"
    found_double = False
    found_single_OH = False
    for nbr in oxy_atoms:
        bond = mol.GetBondBetweenAtoms(Cc_atom.GetIdx(), nbr.GetIdx())
        if bond is None:
            continue
        # Check the bond type.
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            found_double = True
        elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            # For the OH oxygen, ensure that it has a hydrogen attached.
            Hs = [nb for nb in nbr.GetNeighbors() if nb.GetAtomicNum() == 1]
            if len(Hs) >= 1:
                found_single_OH = True
    if not (found_double and found_single_OH):
        return False, "The carboxylic acid group does not appear to be protonated/free (expected C(=O)[OH])"
    
    # --- Check α-carbon chirality (CIP code) ---
    if not Ca_atom.HasProp("_CIPCode"):
        return False, "Alpha-carbon lacks a CIP code; cannot determine configuration"
    cip = Ca_atom.GetProp("_CIPCode")
    if cip != "S":
        return False, f"Alpha-amino acid backbone found but alpha-carbon CIP code is '{cip}', not 'S'"
    
    return True, "Found a free alpha-amino acid backbone with L (S CIP) configuration at the alpha-carbon"

# Example usage:
if __name__ == '__main__':
    # Sample tests from the provided list
    test_smiles = [
        "N[C@@H](CC(=O)c1cc(O)ccc1N)C(O)=O",  # 5-hydroxy-L-kynurenine
        "N[C@@H](CC(=C)C(N)=O)C(O)=O",         # 4-methylene-L-glutamine
        "CC[C@H](N)C(O)=O",                   # L-alpha-aminobutyric acid
        "O=C(O)[C@@H](N)CCC1C=CC(N)C=C1",      # Amiclenomycin
        "N[C@@H](CCCCC(O)=O)C(O)=O",           # L-2-aminopimelic acid
        "N[C@@H](CCCCNC(O)=O)C(O)=O",          # N(6)-carboxy-L-lysine
        # You may add additional examples from the list above.
    ]
    for smi in test_smiles:
        res, reason = is_L_alpha_amino_acid(smi)
        print(smi, "->", res, ":", reason)