"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: L-alpha-amino acid
Definition: Any alpha-amino acid having L-configuration at the alpha-carbon, meaning that the free (non-amidated) amino acid backbone has a free (protonated) carboxylic acid group and the α-carbon (chiral center) is assigned the S CIP configuration.
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
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to allow accurate connectivity inspection.
    mol = Chem.AddHs(mol)
    
    # Compute stereochemistry and assign CIP codes.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define two SMARTS patterns that capture an amino acid backbone.
    # Pattern explanation:
    #  • [NX3;!$(N-C(=O))] : an amino nitrogen not amidated.
    #  • [C@H] or [C@@H] with at least one carbon substituent ([#6])
    #  • A carboxyl group given as [C](=O)[O] where the single-bonded O is expected to be protonated.
    pattern1 = Chem.MolFromSmarts("[NX3;H1,H2;!$(N-C(=O))]-[C@H]([#6])-[C](=O)[O]")
    pattern2 = Chem.MolFromSmarts("[NX3;H1,H2;!$(N-C(=O))]-[C@@H]([#6])-[C](=O)[O]")
    
    matches = []
    for label, pattern in (("pattern1", pattern1), ("pattern2", pattern2)):
        submatches = mol.GetSubstructMatches(pattern, useChirality=True)
        # Store the pattern label with the match for clarity.
        for m in submatches:
            matches.append((label, m))
    
    # Remove duplicate matches (if any) based on atom indices.
    unique_matches = []
    seen = set()
    for label, match in matches:
        tup = tuple(match)
        if tup not in seen:
            seen.add(tup)
            unique_matches.append((label, match))
    
    if len(unique_matches) == 0:
        return False, "Alpha-amino acid backbone (free NH2 and COOH) not found"
    if len(unique_matches) > 1:
        return False, f"Found {len(unique_matches)} amino acid backbone motifs; likely a peptide or multiple motifs"
    
    # Only one backbone motif.
    label, match = unique_matches[0]
    # By convention, our SMARTS is arranged as follows:
    # match[0]: amino nitrogen, match[1]: alpha-carbon, match[2]: carboxyl carbon.
    N_idx, Ca_idx, Cc_idx = match[0], match[1], match[2]
    N_atom = mol.GetAtomWithIdx(N_idx)
    Ca_atom = mol.GetAtomWithIdx(Ca_idx)
    Cc_atom = mol.GetAtomWithIdx(Cc_idx)
    
    # 1. Check that the amino nitrogen is free (should have only one heavy neighbor, the alpha-carbon).
    heavy_neighbors_N = [nbr for nbr in N_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    if len(heavy_neighbors_N) != 1 or heavy_neighbors_N[0].GetIdx() != Ca_idx:
        return False, "Amino nitrogen appears modified or is not free (unexpected heavy atom connectivity)"
    
    # 2. Check the carboxylic acid group.
    # Instead of enforcing exactly 2 oxygen neighbors, we relax to check that:
    # - At least one oxygen is double-bonded to the carboxyl carbon (carbonyl oxygen)
    # - At least one oxygen is single-bonded to the carboxyl carbon and carries at least one hydrogen (hydroxyl oxygen)
    oxy_neighbors = [nbr for nbr in Cc_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(oxy_neighbors) < 2:
        return False, "Carboxyl carbon has fewer than two oxygen neighbors"
    
    found_double = False
    found_protonated_single = False
    for o in oxy_neighbors:
        bond = mol.GetBondBetweenAtoms(Cc_atom.GetIdx(), o.GetIdx())
        if bond is None:
            continue
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            found_double = True
        elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            # Check if the oxygen has at least one hydrogen attached.
            if any(neigh.GetAtomicNum() == 1 for neigh in o.GetNeighbors()):
                found_protonated_single = True
    if not found_double:
        return False, "Carboxylic acid group does not appear to have a carbonyl oxygen (C=O)"
    if not found_protonated_single:
        return False, "Carboxylic acid group does not appear free/protonated (expected -OH)"
    
    # 3. Check chirality of the alpha-carbon.
    if not Ca_atom.HasProp("_CIPCode"):
        return False, "Alpha-carbon lacks a CIP code; cannot determine configuration"
    cip = Ca_atom.GetProp("_CIPCode")
    if cip != "S":
        return False, f"Alpha-carbon CIP code is '{cip}', not 'S' as required for L configuration"
    
    return True, "Found a free L-alpha-amino acid (free NH2 and COOH) with L (S) configuration"

# Example usage (testing a few provided SMILES):
if __name__ == '__main__':
    test_examples = [
        ("5-hydroxy-L-kynurenine", "N[C@@H](CC(=O)c1cc(O)ccc1N)C(O)=O"),
        ("4-methylene-L-glutamine", "N[C@@H](CC(=C)C(N)=O)C(O)=O"),
        ("L-alpha-aminobutyric acid", "CC[C@H](N)C(O)=O"),
        ("Amiclenomycin", "O=C(O)[C@@H](N)CCC1C=CC(N)C=C1"),
        ("L-2-aminopimelic acid", "N[C@@H](CCCCC(O)=O)C(O)=O"),
        ("N(6)-carboxy-L-lysine", "N[C@@H](CCCCNC(O)=O)C(O)=O"),
        ("N(4)-methyl-L-asparagine", "CNC(=O)C[C@H](N)C(O)=O"),
        ("L-serine O-sulfate", "N[C@@H](COS(O)(=O)=O)C(O)=O"),
        ("L-glutamine", "N[C@@H](CCC(N)=O)C(O)=O"),
        ("N(pros)-methyl-L-histidine", "Cn1cncc1C[C@H](N)C(O)=O")
    ]
    for name, smi in test_examples:
        result, reason = is_L_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {result}\nReason: {reason}\n{'-'*65}")