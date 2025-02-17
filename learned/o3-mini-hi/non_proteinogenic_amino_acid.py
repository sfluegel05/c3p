"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
#!/usr/bin/env python
"""
Classifies: Non-proteinogenic amino acid
Definition: Any amino-acid that is not naturally encoded in the genetic code of any organism.

The function is_non_proteinogenic_amino_acid takes a SMILES string as input and returns:
    (bool, str) -> (True, reason) if the molecule appears to be a free amino acid 
                     that is NOT one of the canonical proteinogenic amino acids,
                  (False, reason) otherwise.

The improved algorithm applies these heuristics:
  1. Parse the SMILES and add explicit hydrogens.
  2. Reject molecules that appear to be peptides by the presence of an internal amide bond
     (using a substructure search with SMARTS "C(=O)N[C]").
  3. Identify at least one free carboxyl group using a SMARTS pattern for a carboxyl moiety 
     ([CX3](=O)[OX1,OX2]). We then check that the oxygen in the –OH has at least one hydrogen,
     ensuring that it is not an amide.
  4. Scan for a free amino group (nitrogen atom) that carries at least one hydrogen and is not 
     directly involved in an amide bond (i.e. not bonded to a carbon that bears a double-bonded oxygen).
  5. Finally, compute the InChIKey and reject molecules that match one of the 20 canonical amino acids.
     
Note: This heuristic will not be perfect. Some non–proteinogenic amino acids with unusual connectivity 
(e.g. those lacking an obvious nearby “α-carbon”) might be missed, but this approach avoids falsely 
rejecting candidates solely because the free amino group is not located immediately next to the carboxyl group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# Pre-compute InChIKeys for the 20 canonical proteinogenic amino acids.
_PROTEINOGENIC_SMILES = [
    "NCC(=O)O",                         # Glycine
    "CC(N)C(=O)O",                       # Alanine
    "CC(C)[C@H](N)C(=O)O",                # Valine
    "CCC[C@H](N)C(=O)O",                  # Leucine
    "CC[C@H](C)[C@H](N)C(=O)O",           # Isoleucine
    "C1CC[C@@H](NC1)C(=O)O",              # Proline
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
_PROTEINOGENIC_INCHIKEYS = set()
for smi in _PROTEINOGENIC_SMILES:
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        ik = Chem.MolToInchiKey(mol)
        _PROTEINOGENIC_INCHIKEYS.add(ik)

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines whether the SMILES string represents a free amino acid that is
    non-proteinogenic (i.e. not one of the 20 canonical amino acids).

    The function checks:
      - That the molecule parses and explicit hydrogens are added.
      - That the molecule does not display internal peptide connectivity (an amide
        bond between two amino acid fragments).
      - That the molecule contains at least one free carboxyl group, defined as a
        carbon atom ([CX3]) double-bonded to an oxygen and single-bonded to an oxygen
        atom which carries at least one hydrogen (i.e. an –OH group).
      - That the molecule contains at least one free amino group. A free amino group is
        defined as a nitrogen (atomic number 7) with one or more hydrogens that is not
        involved in an amide bond (i.e. not neighboring a carbonyl group).
      - That the molecule’s InChIKey is not found among the 20 canonical proteinogenic amino acids.
      
    Args:
        smiles (str): SMILES string.

    Returns:
        (bool, str): (True, reason) if classified as a non-proteinogenic amino acid;
                     (False, reason) otherwise.
    """
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens.
    mol = Chem.AddHs(mol)

    # Check for an internal peptide bond. The SMARTS "C(=O)N[C]" looks for a carbonyl
    # connected to a nitrogen that is in turn connected to another carbon – a typical motif
    # in peptides (i.e. between amino acid units).
    peptide_bb = Chem.MolFromSmarts("C(=O)N[C]")
    if mol.HasSubstructMatch(peptide_bb):
        return False, "Molecule appears to be part of a peptide (has internal amide bonds)"
    
    # Look for at least one free carboxyl group.
    # The pattern identifies a carbonyl where the oxygen attached is in an –OH group.
    carboxyl_smarts = Chem.MolFromSmarts("[CX3](=O)[OX1,OX2]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_smarts)
    free_carboxyl_found = False
    for match in carboxyl_matches:
        # For each match, check that at least one oxygen in the pattern (position 1) has a hydrogen.
        oxy_atom = mol.GetAtomWithIdx(match[1])
        if oxy_atom.GetTotalNumHs() >= 1:
            free_carboxyl_found = True
            break
    if not free_carboxyl_found:
        return False, "No free carboxyl group found"
    
    # Look for at least one free amino group.
    # We iterate over nitrogen atoms and require that they have at least one attached hydrogen.
    # We then exclude those nitrogen atoms that are likely involved in an amide bond (i.e.
    # bonded to a carbon that itself is bonded to a double-bonded oxygen).
    free_amino_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        # Check that the nitrogen has at least one hydrogen attached.
        if atom.GetTotalNumHs() < 1:
            continue
        is_amide = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # candidate carbon neighbor
                # Look for a bonded oxygen with a double bond.
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetAtomicNum() == 8 and nbr2.GetSymbol() == "O":
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                        # bond.GetBondType() == Chem.rdchem.BondType.DOUBLE is also acceptable.
                        if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                            is_amide = True
                            break
                if is_amide:
                    break
        if not is_amide:
            free_amino_found = True
            break
    if not free_amino_found:
        return False, "No free amino group found"
    
    # Compute the InChIKey of the hydrogen-added molecule.
    try:
        inchi_key = Chem.MolToInchiKey(mol)
    except Exception as e:
        return False, f"Failed to generate InChIKey: {e}"
    
    if inchi_key in _PROTEINOGENIC_INCHIKEYS:
        return False, "Matches a canonical proteinogenic amino acid"
    
    return True, "Contains a free carboxyl group and a free amino group without internal peptide bonds"

# Example usage (uncomment to test):
# test_smiles = [
#     "[C@H]1(C(N(C1)[C@H](C=2C=CC(=CC2)O)C(=O)O)=O)NC([C@@H](C3=CC=C(C=C3)OCC[C@H](C(O)=O)N)=O",  # nocardicin C - should return False (peptide)
#     "NCCCC[C@H](N)CC(O)=O",  # (3S)-3,7-diaminoheptanoic acid - should return True
#     "NCC(=O)O",             # Glycine - proteinogenic so should return False
#     "CCN(CC)CC(O)=O",        # Example: N,N-diethylglycine, non-proteinogenic, should return True
# ]
# for smi in test_smiles:
#     result, reason = is_non_proteinogenic_amino_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")